# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import re
import os
import copy as cp
from .._utils import run_command
from q2_types_genomics.feature_data import MAGSequencesDirFmt
from q2_types_genomics.genome_data import (
    LociDirectoryFormat, GenesDirectoryFormat, ProteinsDirectoryFormat,
)


def predict_genes_prodigal(
        mags: MAGSequencesDirFmt,
        translation_table_number: int = 11,
) -> (LociDirectoryFormat, GenesDirectoryFormat, ProteinsDirectoryFormat):

    # Instantiate output directories
    loci = LociDirectoryFormat()
    genes = GenesDirectoryFormat()
    proteins = ProteinsDirectoryFormat()

    # Get paths to directories (they already exist as tmp dirs)
    loci_dir = str(loci)
    genes_dir = str(genes)
    proteins_dir = str(proteins)

    # Define base command
    base_cmd = [
        "prodigal",
        "-g", translation_table_number,
        "-f", "gff"
    ]

    # Pattern for input files
    pattern = re.compile(
        r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-4[0-9a-fA-F]{3}-"
        r"[89abAB][0-9a-fA-F]{3}-[0-9a-fA-F]{12}\.(fa|fasta)$"
    )

    # For every fasta file in dir_with_MAGs call prodigal and write
    # outputs corresponding directories.
    for fasta_file in glob.glob(f'{mags.path}/{pattern}'):
        # Get the filename from the file path
        filename = os.path.basename(fasta_file)
        file_id = os.path.splitext(filename)[0]

        # Build paths to outputs TODO: verify regex
        o = os.path.join(loci_dir, f"{file_id}_loci.gff")
        a = os.path.join(proteins_dir, f"{file_id}_proteins.fasta")
        d = os.path.join(genes_dir, f"{file_id}_genes.fasta")

        # Adjust command and run
        cmd = cp.deepcopy(base_cmd)
        cmd.extend([
            "-o", o,
            "-a", a,
            "-d", d
        ])
        run_command(cmd)

    # Return output directories
    return loci, genes, proteins
