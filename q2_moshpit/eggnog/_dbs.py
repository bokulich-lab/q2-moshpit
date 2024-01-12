# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from q2_types.feature_data import ProteinSequencesDirectoryFormat
import shutil
from q2_types_genomics.reference_db import (
    EggnogRefDirFmt, DiamondDatabaseDirFmt, NCBITaxonomyDirFmt,
    EggnogProteinSequencesDirFmt
)
from .._utils import run_command, _process_common_input_params, colorify
from ._utils import _parse_build_diamond_db_params


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


def build_custom_diamond_db(
        seqs: ProteinSequencesDirectoryFormat,
        taxonomy: NCBITaxonomyDirFmt = None,
        threads: int = 1,
        file_buffer_size: int = 67108864,
        ignore_warnings: bool = False,
        no_parse_seqids: bool = False
        ) -> DiamondDatabaseDirFmt:
    '''
    Builds diamond database from protein reference database file in FASTA
    format.
    '''
    # Process input parameters
    kwargs = {}
    for key, value in locals().items():
        if key not in ["seqs", "taxonomy", "kwargs"]:
            kwargs[key] = value

    # Add paths to taxonomy data if provided
    if taxonomy is not None:
        kwargs["taxonmap"] = os.path.join(
            str(taxonomy), "prot.accession2taxid.gz"
        )
        kwargs["taxonnodes"] = os.path.join(str(taxonomy), "nodes.dmp")
        kwargs["taxonnames"] = os.path.join(str(taxonomy), "names.dmp")

    # Filter out all kwargs that are falsy (except 0 and 0.0)
    parsed_args = _process_common_input_params(
        processing_func=_parse_build_diamond_db_params, params=kwargs
    )

    # Instantiate output object
    diamond_db = DiamondDatabaseDirFmt()

    # Run diamond makedb
    cmd = [
        "diamond", "makedb",
        "--verbose", "--log",
        "--in", f"{os.path.join(str(seqs), 'protein-sequences.fasta')}",
        "--db", f"{os.path.join(str(diamond_db), 'ref_db.dmnd')}",
        *parsed_args
    ]
    run_command(cmd)

    # Return output artifact
    return diamond_db


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
    print(colorify("Starting download..."))
    run_command(
        cmd=[
            "wget", "-e", "robots=off", "-O", f"{path_out}",
            "http://eggnogdb.embl.de/download/emapperdb-5.0.2/"
            "eggnog_proteins.dmnd.gz"
        ]
    )

    # Decompressing file
    print(colorify(
            "Download completed.\n"
            "Decompressing file..."
    ))
    run_command(
        cmd=["gunzip", f"{path_out}"]
    )

    # Let user know that the process is done.
    # The actual copying will be taken care of by qiime behind the
    # scenes.
    print(colorify(
        "Decompression completed. \n"
        "Copying file from temporary directory to final location "
        "(this will take a few minutes)..."
    ))

    # Return object
    return diamond_db


def fetch_eggnog_proteins() -> EggnogProteinSequencesDirFmt:
    """
    Downloads eggnog proteome database.
    This script downloads 2 files (e5.proteomes.faa and e5.taxid_info.tsv)
    and creates and artifact with them. At least 18 GB of storage space is
    required to run this action.
    """
    # Initialize output objects
    eggnog_fa = EggnogProteinSequencesDirFmt()
    fasta_file = os.path.join(str(eggnog_fa), "e5.proteomes.faa")
    taxonomy_file = os.path.join(str(eggnog_fa), "e5.taxid_info.tsv")

    # Download fasta file
    print(colorify("Downloading fasta file..."))
    run_command(
        cmd=[
            "wget", "-e", "robots=off", "-O", f"{fasta_file}",
            "http://eggnog5.embl.de/download/eggnog_5.0/e5.proteomes.faa"
        ]
    )

    # Download taxonomy file
    print(colorify(
        "Download completed.\n"
        "Downloading taxonomy file..."
    ))
    run_command(
        cmd=[
            "wget", "-e", "robots=off", "-O", f"{taxonomy_file}",
            "http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv"
        ]
    )

    # Let user know that the process is done.
    # The actual copying will be taken care of by qiime behind the
    # scenes.
    print(colorify(
        "Download completed. \n"
        "Copying files from temporary directory to final location "
        "(this will take a few minutes)..."
    ))

    return eggnog_fa


def build_eggnog_diamond_db(
        eggnog_proteins: EggnogProteinSequencesDirFmt,
        taxon: str
) -> DiamondDatabaseDirFmt:
    """
    Creates an DIAMOND database which contains the protein
    sequences that belong to the specified taxon.
    """

    # Initialize output objects
    diamond_db = DiamondDatabaseDirFmt()

    # Define command.
    cmd = [
        "create_dbs.py",
        "--data_dir", str(eggnog_proteins),
        "--taxids", taxon,
        "--dbname", "ref_db"
    ]
    run_command(cmd)

    # The script will create the diamond DB in side the directory of
    # eggnog_proteins object, so we need to move it to diamond_db
    source_path = os.path.join(str(eggnog_proteins), "ref_db.dmnd")
    destination_path = os.path.join(str(diamond_db), "ref_db.dmnd")
    shutil.move(source_path, destination_path)

    # Return objects
    return diamond_db
