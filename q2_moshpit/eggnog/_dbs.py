# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import shutil
import pandas as pd
from qiime2.core.exceptions import ValidationError
from q2_types.feature_data import ProteinSequencesDirectoryFormat
from q2_types.reference_db import (
    EggnogRefDirFmt, DiamondDatabaseDirFmt, NCBITaxonomyDirFmt,
    EggnogProteinSequencesDirFmt, HmmerDirFmt
)
from q2_moshpit._utils import (
    run_command, _process_common_input_params, colorify,
    _calculate_md5_from_file
)
from q2_moshpit.eggnog._utils import (
    _parse_build_diamond_db_params, _download_and_build_hmm_db,
    _download_fastas_into_hmmer_db, _try_wget
)
import tempfile


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
        taxon: int
) -> DiamondDatabaseDirFmt:
    """
    Creates a DIAMOND database which contains the protein
    sequences that belong to the specified taxon.
    """
    # Validate taxon ID
    _validate_taxon_id(eggnog_proteins, taxon)

    # Initialize output objects
    diamond_db = DiamondDatabaseDirFmt()

    # Define command.
    cmd = [
        "create_dbs.py",
        "--data_dir", str(eggnog_proteins),
        "--taxids", str(taxon),
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


def _validate_taxon_id(eggnog_proteins, taxon):
    # Validate taxon id number
    # Read in valid taxon ids
    taxid_info = pd.read_csv(
        os.path.join(str(eggnog_proteins), "e5.taxid_info.tsv"),
        sep="\t"
    )

    # Convert them into a set
    tax_ids = set()
    for lineage in taxid_info["Taxid Lineage"]:
        tax_ids.update(
            set(
                lineage.strip().split(",")
            )
        )

    # Check for overlap with provided taxon id
    if not str(taxon) in tax_ids:
        raise ValueError(
            f"'{taxon}' is not valid taxon ID. "
            "To view all valid taxon IDs inspect e5.taxid_info.tsv. "
            "You can download it with this command: "
            "wget "
            "http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv"
        )


def fetch_ncbi_taxonomy() -> NCBITaxonomyDirFmt:
    """
    Script fetches 3 files from the NCBI server and puts them into the folder
    of a NCBITaxonomyDirFmt object.
    """
    ncbi_data = NCBITaxonomyDirFmt()
    zip_path = os.path.join(str(ncbi_data), "taxdmp.zip")
    proteins_path = os.path.join(str(ncbi_data), "prot.accession2taxid.gz")

    # Download dump zip file + MD5 file
    print(colorify("Downloading *.dmp files..."))
    run_command(
        cmd=[
            "wget", "-O", f"{zip_path}",
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
        ]
    )
    run_command(
        cmd=[
            "wget", "-O", f"{zip_path}.md5",
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip.md5"
        ]
    )

    _collect_and_compare_md5(f"{zip_path}.md5", zip_path)

    run_command(
        cmd=[
            "unzip", "-j", zip_path, "names.dmp", "nodes.dmp",
            "-d", str(ncbi_data)
        ]
    )

    os.remove(zip_path)

    # Download proteins + MD5 file
    print(colorify("Downloading proteins file (~8 GB)..."))
    run_command(
        cmd=[
            "wget", "-O", f"{proteins_path}",
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
            "prot.accession2taxid.gz"
        ]
    )
    run_command(
        cmd=[
            "wget", "-O", f"{proteins_path}.md5",
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
            "prot.accession2taxid.gz.md5"
        ]
    )

    _collect_and_compare_md5(f"{proteins_path}.md5", proteins_path)

    print(colorify(
        "Done! Moving data from temporary directory to final location..."
    ))
    return ncbi_data


def _collect_and_compare_md5(path_to_md5: str, path_to_file: str):
    # Read in hash from md5 file
    with open(path_to_md5, 'r') as f:
        expected_hash = f.readline().strip().split(maxsplit=1)[0]

    # Calculate hash from file
    observed_hash = _calculate_md5_from_file(path_to_file)

    if observed_hash != expected_hash:
        raise ValidationError(
            "Download error. Data possibly corrupted.\n"
            f"{path_to_file} has an unexpected MD5 hash.\n\n"
            "Expected hash:\n"
            f"{expected_hash}\n\n"
            "Observed hash:\n"
            f"{observed_hash}"
        )

    # If no exception is raised, remove md5 file
    os.remove(path_to_md5)


def fetch_eggnog_hmmer_db(taxon_id: int) -> HmmerDirFmt:
    # Validate taxon ID
    with tempfile.TemporaryDirectory() as tmp:
        print(colorify(
                "Validating taxon ID: \n"
                "Downloading taxonomy file..."
        ))
        _try_wget(
            f"{tmp}/e5.taxid_info.tsv",
            "http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv",
            "Error during taxon-info-file download"
        )
        _validate_taxon_id(tmp, taxon_id)

    # Download HMMER data base
    print(colorify(
        "Valid taxon ID. \n"
        "Proceeding with HMMER database download and build..."
    ))
    hmmer_db = _download_and_build_hmm_db(taxon_id)
    print(colorify(
        "HMM database built successfully. \n"
        "Proceeding with FASTA files download and processing..."
    ))

    # Downland fasta sequences
    _download_fastas_into_hmmer_db(hmmer_db, taxon_id)
    print(colorify(
        "FASTA files processed successfully. \n"
        "Moving data from temporary to final location..."
    ))

    return hmmer_db
