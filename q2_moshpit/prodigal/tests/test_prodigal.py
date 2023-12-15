# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from q2_moshpit.prodigal.prodigal import predict_genes_prodigal
from qiime2.plugin.testing import TestPluginBase
from q2_types_genomics.feature_data import MAGSequencesDirFmt
from unittest.mock import patch, call
from q2_types_genomics.genome_data import (
    LociDirectoryFormat, GenesDirectoryFormat, ProteinsDirectoryFormat,
)


class TestBUSCO(TestPluginBase):
    package = "q2_moshpit.prodigal.tests"

    @patch("subprocess.run")
    def test_run_prodigal_1_mag(self, subp_run):
        # Run prodigal with dummy data
        p = self.get_data_path("dir_with_1_mag")
        mags = MAGSequencesDirFmt(path=p, mode="r")
        loci, genes, proteins = predict_genes_prodigal(mags=mags)

        # Check that output is correct type
        self.assertIsInstance(loci, LociDirectoryFormat)
        self.assertIsInstance(genes, GenesDirectoryFormat)
        self.assertIsInstance(proteins, ProteinsDirectoryFormat)

        # Get names of fasta files from test data dir
        fasta_file = [
            os.path.splitext(file)[0] for file in os.listdir(mags.path)
            if file.endswith(".fa") or file.endswith(".fasta")
        ]

        # Unpack list
        fasta_file = fasta_file[0]

        # Assert that patch was called once
        subp_run.assert_called_once_with([
            "prodigal",
            "-g", "11",
            "-f", "gff",
            "-i", os.path.join(mags.path, f"{fasta_file}.fasta"),
            "-o", os.path.join(loci.path, f"{fasta_file}_loci.gff"),
            "-a", os.path.join(proteins.path, f"{fasta_file}_proteins.fasta"),
            "-d", os.path.join(genes.path, f"{fasta_file}_genes.fasta")],
            cwd=None,
            check=True
        )

    @patch("subprocess.run")
    def test_run_prodigal_3_mag(self, subp_run):
        # Run prodigal with dummy data
        p = self.get_data_path("dir_with_3_mag")
        mags = MAGSequencesDirFmt(path=p, mode="r")
        loci, genes, proteins = predict_genes_prodigal(mags=mags)

        # Check that output is correct type
        self.assertIsInstance(loci, LociDirectoryFormat)
        self.assertIsInstance(genes, GenesDirectoryFormat)
        self.assertIsInstance(proteins, ProteinsDirectoryFormat)

        # Get names of fasta files from test data dir
        fasta_files = [
            os.path.splitext(file)[0] for file in os.listdir(mags.path)
            if file.endswith(".fa") or file.endswith(".fasta")
        ]

        # Define calls
        three_calls = [call([
            "prodigal",
            "-g", "11",
            "-f", "gff",
            "-i", os.path.join(mags.path, f"{fasta_file}.fasta"),
            "-o", os.path.join(loci.path, f"{fasta_file}_loci.gff"),
            "-a", os.path.join(proteins.path, f"{fasta_file}_proteins.fasta"),
            "-d", os.path.join(genes.path, f"{fasta_file}_genes.fasta")],
             cwd=None, check=True)
            for fasta_file in fasta_files
        ]

        # Assert that patch was called 3 times
        subp_run.assert_has_calls(three_calls)
