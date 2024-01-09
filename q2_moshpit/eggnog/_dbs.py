# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from .._utils import run_command, colorify
from q2_types_genomics.reference_db import (
    EggnogRefDirFmt, DiamondDatabaseDirFmt
)


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
    # y: Answer yes to all prompts thrown by download_eggnog_data.py
    # D: Do not download the Diamond database
    # data_dir: location where to save downloads
    cmd = [
        "download_eggnog_data.py", "-y", "-D",
        "--data_dir", str(eggnog_db.path)
    ]
    run_command(cmd)

    # Return objects
    return eggnog_db


def fetch_diamond_db() -> DiamondDatabaseDirFmt:
    """
    Downloads diamond reference database using the
    `download_eggnog_data.py` script from eggNOG. Here, this
    script downloads 1 file (8.6 Gb).
    """

    # Initialize output objects
    diamond_db = DiamondDatabaseDirFmt()
    path_out = os.path.join(str(diamond_db), "ref_db.dmnd.gz")

    # Download Diamond DB
    print(
        colorify(
            "Starting download...", "lgreen"
        )
    )
    run_command(
        cmd=[
            "wget", "-e", "robots=off", "-O", f"{path_out}",
            "http://eggnogdb.embl.de/download/emapperdb-5.0.2/"
            "eggnog_proteins.dmnd.gz"
        ]
    )

    # Decompressing file
    print(
        colorify(
            "Download completed.\n"
            "Decompressing file...",
            "lgreen"
        )
    )
    run_command(
        cmd=["gunzip", f"{path_out}"]
    )

    # Let user know that the process is done.
    # The actual copying wil be taken care of by qiime behind the
    # scenes.
    print(
        colorify(
            "Decompression completed. \n"
            "Copying file from temporary directory to final location "
            "(this will take a few minutes)...",
            "lgreen"
        )
    )

    # Return object
    return diamond_db
