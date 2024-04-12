# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import copy as cp
from typing import Union
from .._utils import run_command
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt
from q2_types.genome_data import (
    LociDirectoryFormat, GenesDirectoryFormat, ProteinsDirectoryFormat,
)


def predict_genes_prodigal(
        mags: Union[MAGSequencesDirFmt, MultiMAGSequencesDirFmt],
        translation_table_number: str = "11",
) -> (LociDirectoryFormat, GenesDirectoryFormat, ProteinsDirectoryFormat):

    # Instantiate output directories
    loci = LociDirectoryFormat()
    genes = GenesDirectoryFormat()
    proteins = ProteinsDirectoryFormat()

    # Define base command
    base_cmd = [
        "prodigal",
        "-g", translation_table_number,
        "-f", "gff"
    ]

    def _process_fasta_files(fasta_files: list, prefix: str, input_path: str):
        # For every fasta file call prodigal and write
        # outputs to the corresponding directories.
        for fasta_file in fasta_files:
            # Get the filename from the file path
            file_id = os.path.splitext(fasta_file)[0]

            # Adjust command and run
            cmd = cp.deepcopy(base_cmd)
            cmd.extend([
                "-i", os.path.join(input_path, fasta_file),
                "-o",
                os.path.join(
                    loci.path, f"{prefix}{file_id}_loci.gff"
                ),
                "-a",
                os.path.join(
                    proteins.path, f"{prefix}{file_id}_proteins.fasta"
                ),
                "-d",
                os.path.join(
                    genes.path, f"{prefix}{file_id}_genes.fasta"
                )
            ])
        run_command(cmd)

    if isinstance(mags, MAGSequencesDirFmt):
        # Get paths to fasta files in input dir
        fasta_files = os.listdir(mags.path)
        _process_fasta_files(fasta_files, '', mags.path)

    elif isinstance(mags, MultiMAGSequencesDirFmt):
        # List all directories / samples
        for sample_dir in os.listdir(mags.path):
            if os.path.isdir(os.path.join(mags.path, sample_dir)):
                fasta_files = os.listdir(os.path.join(mags.path, sample_dir))
                _process_fasta_files(
                    fasta_files,
                    f"{sample_dir}_",
                    os.path.join(mags.path, sample_dir)
                )

    # Return output directories
    return loci, genes, proteins
