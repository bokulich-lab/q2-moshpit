# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from q2_types_genomics.reference_db import (
    EggnogRefDirFmt, EggnogSequenceTaxaDirFmt
)
from .._utils import run_command, colorify


def fetch_eggnog_db() -> EggnogRefDirFmt:
    """
    Downloads eggnog reference database using the
    `download_eggnog_data.py` script from eggNOG. Here, this
    script downloads 3 files amounting to 47Gb in total.
    """

    # Initialize output objects
    eggnog_db = EggnogRefDirFmt()

    # Define command.
    # Meaning of flags:
    # y: Answer yest to all prompts thrown by download_eggnog_data.py
    # D: Do not download the Diamond database
    # data_dir: location where to save downloads
    cmd = [
        "download_eggnog_data.py", "-y", "-D",
        "--data_dir", str(eggnog_db.path)
    ]
    run_command(cmd)

    # Return objects
    return eggnog_db


def fetch_eggnog_fasta() -> EggnogSequenceTaxaDirFmt:
    """
    # TODO: Add description
    """
    # Initialize output objects
    eggnog_fa = EggnogSequenceTaxaDirFmt()
    fasta_file = os.path.join(str(eggnog_fa), "e5.proteomes.faa")
    taxonomy_file = os.path.join(str(eggnog_fa), "e5.taxid_info.tsv")

    # Download Diamond DB
    print(
        colorify(
            "Downloading fasta file...", "lgreen"
        )
    )
    run_command(
        cmd=[
            "wget", "-e", "robots=off", "-O", f"{fasta_file}",
            "http://eggnog5.embl.de/download/eggnog_5.0/e5.proteomes.faa"
        ]
    )

    # Decompressing file
    print(
        colorify(
            "Download completed.\n"
            "Downloading taxonomy file...",
            "lgreen"
        )
    )
    run_command(
        cmd=[
            "wget", "-e", "robots=off", "-O", f"{taxonomy_file}",
            "http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv"
        ]
    )

    # Let user know that the process is done.
    # The actual copying wil be taken care of by qiime behind the
    # scenes.
    print(
        colorify(
            "Download completed. \n"
            "Copying files from temporary directory to final location "
            "(this will take a few minutes)...",
            "lgreen"
        )
    )

    return eggnog_fa
