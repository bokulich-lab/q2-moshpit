# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
from copy import deepcopy
from typing import List

import requests
import xmltodict
from q2_types.feature_data import DNAFASTAFormat

from q2_moshpit._utils import _process_common_input_params, run_command
from q2_types_genomics.kraken2 import (
    Kraken2DBDirectoryFormat, BrackenDBDirectoryFormat,
)

from q2_moshpit.kraken2.utils import _process_kraken2_arg


COLLECTIONS = {
    "standard": "standard",
    "viral": "viral",
    "minusb": "minusb",
    "standard8": "standard_08gb",
    "standard16": "standard_16gb",
    "pluspf": "pluspf",
    "pluspf8": "pluspf_08gb",
    "pluspf16": "pluspf_16gb",
    "pluspfp": "pluspfp",
    "pluspfp8": "pluspfp_08gb",
    "pluspfp16": "pluspfp_16gb",
    "eupathdb": "eupathdb48",
}
S3_COLLECTIONS_URL = 'https://genome-idx.s3.amazonaws.com'

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
            "An error was encountered while building the database, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def _find_latest_db(collection: str, response: requests.Response) -> str:
    collection_id = COLLECTIONS[collection]
    pattern = fr'kraken\/k2_{collection_id}_\d{{8}}.tar.gz'

    s3_objects = xmltodict.parse(response.content)['ListBucketResult']['Contents']
    s3_objects = [obj for obj in s3_objects if re.match(pattern, obj['Key'])]
    s3_objects = sorted(s3_objects, key=lambda x: x['LastModified'], reverse=True)
    latest_db = s3_objects[0]['Key']
    return latest_db


def _fetch_db_collection(collection: str, tmp_dir: str):
    err_msg = 'Could not connect to the server. Please check your internet ' \
              'connection and try again. The error was: {}.'
    try:
        response = requests.get(S3_COLLECTIONS_URL)
    except requests.exceptions.ConnectionError as e:
        raise ValueError(err_msg.format(e))

    if response.status_code == 200:
        latest_db = _find_latest_db(collection, response)
    else:
        raise ValueError(
            'Could not fetch the list of available databases. '
            f'Status code was: {response.status_code}. '
            'Please try again later.'
        )

    db_uri = f'{S3_COLLECTIONS_URL}/{latest_db}'
    try:
        response = requests.get(db_uri)
    except requests.exceptions.ConnectionError as e:
        raise ValueError(err_msg.format(e))

    db_path = os.path.join(tmp_dir, os.path.split(db_uri)[-1])
    with open(db_path, "wb") as f:
        f.write(response.content)
    with tarfile.open(db_path, "r:gz") as tar:
        tar.extractall(path=tmp_dir)


def _move_db_files(source: str, destination: str, extension: str = "k2d"):
    files = glob.glob(f"{source}/*.{extension}")
    for file in files:
        new_file = os.path.join(destination, os.path.split(file)[-1])
        shutil.move(file, new_file)


def build_kraken_db(
    seqs: DNAFASTAFormat = None,
    collection: str = None,
    threads: int = 1,
    use_ftp: bool = False,
    no_masking: bool = False,
) -> (Kraken2DBDirectoryFormat, BrackenDBDirectoryFormat):
    kraken2_db = Kraken2DBDirectoryFormat()
    bracken_db = BrackenDBDirectoryFormat()

    with tempfile.TemporaryDirectory() as tmp:
        # Construct the custom-made database
        if seqs:
            # Fetch taxonomy (also needed for custom databases)
            _fetch_taxonomy(db_dir=tmp, threads=threads, use_ftp=use_ftp)

            for seq in seqs:
                _add_seqs_to_library(
                    db_dir=tmp, seqs=seq, no_masking=no_masking
                )

            # Build the database
            _build_database(db_dir=tmp, all_kwargs=locals())

            # Move the database files (*.k2d) to the final location
            _move_db_files(source=tmp, destination=str(kraken2_db.path))

            # TODO: implement Bracken DB part
        elif collection:
            # Find files with the latest version
            _fetch_db_collection(collection=collection, tmp_dir=tmp)

            # Move the Kraken2/Bracken database files to the final location
            _move_db_files(tmp, str(kraken2_db.path), extension="k2d")
            _move_db_files(tmp, str(bracken_db.path), extension="kmer_distrib")

        else:
            raise ValueError(
                'You need to either provide a list of sequences to build the '
                'database from or a valid collection name to be fetched from '
                '"Kraken 2/Bracken Refseq indexes" resource.'
            )

    return kraken2_db, bracken_db
