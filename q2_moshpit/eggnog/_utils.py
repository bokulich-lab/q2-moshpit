# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import re
import os
import gzip
from tqdm import tqdm
import tempfile
import subprocess
from typing import List
from q2_types.reference_db import HmmerDirFmt
from .._utils import run_command, colorify


def _parse_build_diamond_db_params(arg_key, arg_val) -> List[str]:
    """Creates a list with argument and its value to be consumed by
    the `diamond makedb` command.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.
    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """
    # Change "_" in arg_key for "-"
    arg_key = arg_key.replace("_", "-")

    if isinstance(arg_val, bool):
        return [f"--{arg_key}"]
    else:
        return [f"--{arg_key}", str(arg_val)]


def _download_and_build_hmm_db(taxon_id) -> HmmerDirFmt:
    hmmer_db = HmmerDirFmt()
    with tempfile.TemporaryDirectory() as tmp:
        _try_wget(
            f"{tmp}/{taxon_id}_hmms.tar.gz",
            "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/"
            f"{taxon_id}/{taxon_id}_hmms.tar.gz",
            "Error during HMMER database download"
        )

        # Extracting
        print(colorify("Decompressing..."))
        run_command(
            cmd=["tar", "zxf", f"{taxon_id}_hmms.tar.gz"],
            cwd=tmp
        )

        # Merge hmm files + write .hmm.idmap
        print(colorify("Merging hmm files..."))
        hmms_merged_p = f"{str(hmmer_db)}/{taxon_id}.hmm"
        idmap_p = f"{str(hmmer_db)}/{taxon_id}.hmm.idmap"

        # Open output files
        with open(hmms_merged_p, "a") as hmms, open(idmap_p, "a") as idmap:

            # Iterate through all decompressed files
            for root, dirnames, files in os.walk(f"{tmp}/{taxon_id}"):
                for i, file in tqdm(enumerate(files, start=1)):
                    if file.endswith(".hmm"):

                        # process hmm files
                        with open(f"{root}/{file}", "r") as hmm_file:
                            lines = hmm_file.readlines()

                            # Find "NAME" line
                            for j, line in enumerate(lines):
                                if line.startswith("NAME "):
                                    modified_line = re.sub(
                                        r"\.faa\.final_tree(\.fa)?", "", line
                                    )

                                    # write modified content to hmms_merged
                                    lines[j] = modified_line
                                    hmms.writelines(lines)

                                    # get name and write to idmap
                                    id = modified_line.replace("NAME  ", "", 1)
                                    idmap.write(f"{i} {id}")

                                    break

    # prepare an HMM database for faster hmmscan searches
    print(colorify("Preparing HMM database..."))
    run_command(cmd=["hmmpress", hmms_merged_p])
    os.remove(hmms_merged_p)

    return hmmer_db


def _download_fastas_into_hmmer_db(hmmer_db: HmmerDirFmt, taxon_id: int):
    with tempfile.TemporaryDirectory() as tmp:
        _try_wget(
            f"{tmp}/{taxon_id}_raw_algs.tar",
            "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/"
            f"{taxon_id}/{taxon_id}_raw_algs.tar",
            "Error downloading FASTA files"
        )

        # Extracting
        print(colorify("Decompressing..."))
        run_command(
            cmd=["tar", "xf", f"{taxon_id}_raw_algs.tar"],
            cwd=tmp
        )

        files = [
            f"{tmp}/{taxon_id}/{f}"
            for f in os.listdir(f"{tmp}/{taxon_id}")
            if f.endswith(".gz")
        ]

        # Extract, remove '-' and save to hmmer_db location
        print(colorify("Processing FASTA files (this can take a while)... "))
        for fpi in tqdm(files):
            new_name = os.path.basename(fpi).replace(".raw_alg.faa.gz", ".fa")
            fpo = os.path.join(str(hmmer_db), new_name)
            with gzip.open(fpi, "rt") as f_in, open(fpo, "w") as f_out:
                content = f_in.read()
                content = content.replace("-", "")
                f_out.write(content)


def _try_wget(output_file: str, url: str, exception_msg: str):
    try:
        run_command(cmd=["wget", "-O", output_file, url])
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"{exception_msg}: {e.returncode}"
        )
