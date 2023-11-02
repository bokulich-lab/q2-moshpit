# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_moshpit.prodigal.prodigal import predict_genes_prodigal
from qiime2.plugin.testing import TestPluginBase
from q2_types_genomics.feature_data import MAGSequencesDirFmt
from unittest.mock import patch, ANY, call
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

        # Assert that patch was called once
        subp_run.assert_called_once_with(
            [
                "prodigal",
                "-g", "11",
                "-f", "gff",
                "-i", ANY,
                "-o", ANY,
                "-a", ANY,
                "-d", ANY
            ],
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

        # Define calls
        three_calls = [
            call(
                [
                    "prodigal",
                    "-g", "11",
                    "-f", "gff",
                    "-i", ANY,
                    "-o", ANY,
                    "-a", ANY,
                    "-d", ANY
                ],
                check=True
            ) for i in range(3)
        ]

        # Assert that patch was called 3 times
        subp_run.assert_has_calls(three_calls)

    # Test weather exception is raised when the empty folder is used for input
    @patch("subprocess.run")
    def test_no_mags_in_input(self, subp_run):
        # Run busco and save paths to run summaries
        p = self.get_data_path("dir_with_0_mag")
        mags = MAGSequencesDirFmt(path=p, mode="r")

        with self.assertRaises(FileNotFoundError):
            predict_genes_prodigal(mags=mags)

        #  Assert that prodigal was not run
        subp_run.assert_not_called()
