# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
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

    # Pattern of fasta files in input dir
    pattern = re.compile(
        r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-4[0-9a-fA-F]{3}-"
        r"[89abAB][0-9a-fA-F]{3}-[0-9a-fA-F]{12}\.(fa|fasta)$"
    )

    # Look for fasta files in input dir
    fasta_files = [
        file for file in os.listdir(mags.path) if pattern.match(file)
    ]

    # Raise exception if no fasta files are found
    if len(fasta_files) == 0:
        raise FileNotFoundError(f"No fasta files found in {mags.path}")

    # Define base command
    base_cmd = [
        "prodigal",
        "-g", str(translation_table_number),
        "-f", "gff"
    ]

    # For every fasta file in mags.path call prodigal and write
    # outputs corresponding directories.
    for fasta_file in fasta_files:
        # Get the filename from the file path
        file_id = os.path.splitext(fasta_file)[0]

        # Build paths to outputs
        i = os.path.join(mags.path, fasta_file)
        o = os.path.join(loci.path, f"{file_id}_loci.gff")
        a = os.path.join(proteins.path, f"{file_id}_proteins.fasta")
        d = os.path.join(genes.path, f"{file_id}_genes.fasta")

        # Adjust command and run
        cmd = cp.deepcopy(base_cmd)
        cmd.extend([
            "-i", i,
            "-o", o,
            "-a", a,
            "-d", d
        ])
        run_command(cmd)

    # Return output directories
    return loci, genes, proteins