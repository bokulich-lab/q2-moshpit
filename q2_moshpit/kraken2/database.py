# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import shutil
import subprocess
import tempfile
from contextlib import ExitStack
from copy import deepcopy
from typing import List

from q2_types.feature_data import DNAFASTAFormat

from q2_moshpit._utils import _process_common_input_params, run_command
from q2_types_genomics.kraken2 import (
    Kraken2DBDirectoryFormat,
)

from q2_moshpit.kraken2.utils import _process_kraken2_arg


def _build_standard_db(db_dir: str, all_kwargs: dict):
    kwargs = {
        k: v for k, v in all_kwargs.items()
        if k in ["threads", "no_masking", "use_ftp", "fast_build"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    if all_kwargs["max_db_size"] > 0:
        common_args.extend(
            ["--max-db-size", str(all_kwargs["max_db_size"])]
        )
    cmd = [
        "kraken2-build", "--standard", "--db", db_dir, *common_args
    ]
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while building the standard "
            f"library, (return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _fetch_taxonomy(db_dir: str, threads: int, use_ftp: bool):
    cmd = [
        "kraken2-build", "--download-taxonomy",
        "--threads", str(threads), "--db", str(db_dir),
    ]
    cmd.append("--use-ftp") if use_ftp else False
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while downloading taxonomy, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _fetch_libraries(db_dir: str, libraries: List[str], all_kwargs: dict):
    kwargs = {
        k: v for k, v in all_kwargs.items()
        if k in ["threads", "no_masking", "use_ftp"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    common_args.extend(["--db", db_dir])
    base_cmd = [
        "kraken2-build", "--download-library",
    ]
    fetch_behaviour = all_kwargs.get("library_exists", "refetch")
    for library in libraries:
        if fetch_behaviour == "skip":
            lib_path = os.path.join(db_dir, "library", library)
            if os.path.exists(lib_path) and len(os.listdir(lib_path)) > 0:
                print(
                    f"Skipping download of the '{library}' library, "
                    f"already exists."
                )
                continue
        try:
            cmd = deepcopy(base_cmd)
            cmd.extend([library, *common_args])
            run_command(cmd=cmd, verbose=True)
        except subprocess.CalledProcessError as e:
            raise Exception(
                f"An error was encountered while downloading the "
                f"'{library}' library, (return code {e.returncode}), "
                "please inspect stdout and stderr to learn more."
            )


def _add_seqs_to_library(db_dir: str, seqs: DNAFASTAFormat, no_masking: bool):
    cmd = [
        "kraken2-build", "--add-to-library",
        str(seqs.path), "--db", db_dir
    ]
    cmd.append("--no-mask") if no_masking else False
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while adding sequences to the "
            f"library, (return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _build_database(db_dir: str, all_kwargs: dict):
    kwargs = {
        k: v for k, v in all_kwargs.items()
        if k in ["threads", "minimizer_len", "minimizer_spaces",
                 "load_factor", "fast_build", "kmer_len"]
    }
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    if all_kwargs["max_db_size"] > 0:
        common_args.extend(
            ["--max-db-size", str(all_kwargs["max_db_size"])]
        )
    cmd = [
        "kraken2-build", "--build", "--db", db_dir, *common_args
    ]
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while building the standard  "
            f"library, (return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _move_db_files(source: str, destination: str):
    files = glob.glob(f"{source}/*.k2d")
    for file in files:
        new_file = os.path.join(destination, os.path.split(file)[-1])
        shutil.move(file, new_file)


def build_kraken_db(
    seqs: DNAFASTAFormat = None,
    standard: bool = False,
    library_path: str = None,
    libraries: List[str] = None,
    library_exists: str = 'skip',
    threads: int = 1,
    kmer_len: int = 35,
    minimizer_len: int = 31,
    minimizer_spaces: int = 7,
    no_masking: bool = False,
    max_db_size: int = 0,
    use_ftp: bool = False,
    load_factor: float = 0.7,
    fast_build: bool = False,
) -> Kraken2DBDirectoryFormat:
    db = Kraken2DBDirectoryFormat()
    db_dir = str(db.path)

    if standard and libraries:
        raise ValueError(
            'Standard Kraken2 database was requested but some libraries '
            'were also provided. Please provide either only the "standard" '
            'option or a list of "libraries" to be fetched for the database.'
        )

    with ExitStack() as stack:
        if not library_path:
            temp = tempfile.TemporaryDirectory()
            temp_dir = temp.name
            stack.enter_context(temp)
        else:
            os.makedirs(library_path, exist_ok=True)
            temp_dir = library_path

        # If requested, build the standard Kraken2 database
        if standard:
            _build_standard_db(db_dir=temp_dir, all_kwargs=locals())
            _move_db_files(source=temp_dir, destination=db_dir)
            return db

        # Fetch taxonomy (required regardless of the build source)
        _fetch_taxonomy(db_dir=temp_dir, threads=threads, use_ftp=use_ftp)

        # If requested, download all the libraries
        if libraries:
            _fetch_libraries(
                db_dir=temp_dir, libraries=libraries, all_kwargs=locals()
            )

        # If provided, add the provided sequences to the database
        if seqs:
            for seq in seqs:
                _add_seqs_to_library(
                    db_dir=temp_dir, seqs=seq, no_masking=no_masking
                )

        # Finally, build the actual database
        _build_database(db_dir=temp_dir, all_kwargs=locals())

        # Move the database files (*.k2d) to the final location
        _move_db_files(source=temp_dir, destination=db_dir)

    return db
