# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from q2_types_genomics.reference_db import EggnogRefDirFmt
from q2_types.feature_data import ProteinSequencesDirectoryFormat
from q2_types_genomics.reference_db import (
    DiamondDatabaseDirFmt, NCBITaxonomyDirFmt
)
from .._utils import run_command, _process_common_input_params
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


def build_custom_diamond_db(
        seqs: ProteinSequencesDirectoryFormat,
        taxonomy: NCBITaxonomyDirFmt = None,
        threads: int = 1,
        log: bool = False,
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
        "--verbose",
        "--in", f"{os.path.join(str(seqs), 'protein-sequences.fasta')}",
        "--db", f"{os.path.join(str(diamond_db), 'ref_db.dmnd')}",
        *parsed_args
    ]
    run_command(cmd)

    # Return output artifact
    return diamond_db
