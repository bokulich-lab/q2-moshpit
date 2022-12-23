# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import unittest

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
)
from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
)

from q2_types_genomics.per_sample_data._format import MultiMAGSequencesDirFmt
from unittest.mock import patch, ANY, call

from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.kraken2.kraken2 import (
    _classify_kraken,
    _get_seq_paths,
    _construct_output_paths,
)


class TestKraken2(TestPluginBase):
    package = "q2_moshpit.kraken2.tests"

    def setUp(self):
        super().setUp()

    def test_get_seq_paths_reads_single(self):
        manifest = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path("single-end"), mode="r"
        ).manifest.view(pd.DataFrame)

        expected = [
            (f"sample{i}", f"sample{i}", f"reads{i}_R1.fastq.gz")
            for i in range(1, 3)
        ]
        for ((index, row), (exp_bin, exp_sample, exp_fp)) in zip(
            manifest.iterrows(), expected
        ):
            _sample, _bin, fn = _get_seq_paths(index, row, manifest.columns)
            self.assertEqual(_sample, exp_sample)
            self.assertEqual(_bin, exp_bin)
            self.assertTrue(fn[0].endswith(exp_fp))

    def test_get_seq_paths_reads_paired(self):
        manifest = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path("paired-end"), mode="r"
        ).manifest.view(pd.DataFrame)

        expected = [
            (
                f"sample{i}",
                f"sample{i}",
                [f"reads{i}_R1.fastq.gz", f"reads{i}_R2.fastq.gz"],
            )
            for i in range(1, 3)
        ]
        for ((index, row), (exp_bin, exp_sample, exp_fp)) in zip(
            manifest.iterrows(), expected
        ):
            _sample, _bin, fn = _get_seq_paths(index, row, manifest.columns)
            self.assertEqual(_sample, exp_sample)
            self.assertEqual(_bin, exp_bin)
            self.assertTrue(fn[0].endswith(exp_fp[0]))
            self.assertTrue(fn[1].endswith(exp_fp[1]))

    def test_get_seq_paths_mags(self):
        manifest = MultiMAGSequencesDirFmt(
            self.get_data_path("mags"), "r"
        ).manifest.view(pd.DataFrame)

        expected = [
            ("samp1", "bin1", "samp1/bin1.fa"),
            ("samp1", "bin2", "samp1/bin2.fa"),
            ("samp2", "bin1", "samp2/bin1.fa"),
        ]
        for ((index, row), (exp_sample, exp_bin, exp_fp)) in zip(
            manifest.iterrows(), expected
        ):
            _sample, _bin, fn = _get_seq_paths(index, row, manifest.columns)
            self.assertEqual(_sample, exp_sample)
            self.assertEqual(_bin, exp_bin)
            self.assertTrue(fn[0].endswith(exp_fp))

    @patch("os.makedirs")
    def test_construct_output_paths(self, p1):
        _sample, _bin = "sample1", "bin1"
        reports_dir = Kraken2ReportDirectoryFormat()
        outputs_dir = Kraken2OutputDirectoryFormat()

        exp_rep_fp = os.path.join(
            reports_dir.path, _sample, f"{_bin}.report.txt"
        )
        exp_out_fp = os.path.join(
            outputs_dir.path, _sample, f"{_bin}.output.txt"
        )
        obs_out_fp, obs_rep_fp = _construct_output_paths(
            _sample, _bin, outputs_dir, reports_dir
        )
        self.assertEqual(obs_rep_fp, exp_rep_fp)
        self.assertEqual(obs_out_fp, exp_out_fp)
        p1.assert_has_calls(
            [
                call(os.path.split(exp_rep_fp)[0], exist_ok=True),
                call(os.path.split(exp_out_fp)[0], exist_ok=True),
            ]
        )

    @patch("q2_moshpit.kraken2.Kraken2OutputDirectoryFormat")
    @patch("q2_moshpit.kraken2.Kraken2ReportDirectoryFormat")
    @patch("q2_moshpit.kraken2._get_seq_paths")
    @patch("q2_moshpit.kraken2._construct_output_paths")
    @patch("q2_moshpit.kraken2.run_command")
    def test_classify_kraken(self, p1, p2, p3, p4, p5):
        manifest = MultiMAGSequencesDirFmt(
            self.get_data_path("mags"), "r"
        ).manifest.view(pd.DataFrame)
        common_args = ["--db", "/some/where/db", "--quick"]

        fake_report_dir = Kraken2ReportDirectoryFormat()
        fake_output_dir = Kraken2OutputDirectoryFormat()
        exp_out_fps = [
            os.path.join(fake_output_dir.path, "samp1", "bin1.output.txt"),
            os.path.join(fake_output_dir.path, "samp1", "bin2.output.txt"),
            os.path.join(fake_output_dir.path, "samp2", "bin1.output.txt"),
        ]
        exp_rep_fps = [
            os.path.join(fake_report_dir.path, "samp1", "bin1.report.txt"),
            os.path.join(fake_report_dir.path, "samp1", "bin2.report.txt"),
            os.path.join(fake_report_dir.path, "samp2", "bin1.report.txt"),
        ]

        p2.side_effect = list(zip(exp_out_fps, exp_rep_fps))
        p3.side_effect = [
            ("samp1", "bin1", ["samp1/bin1.fa"]),
            ("samp1", "bin2", ["samp1/bin2.fa"]),
            ("samp2", "bin1", ["samp2/bin1.fa"]),
        ]
        p4.return_value = fake_report_dir
        p5.return_value = fake_output_dir

        # run kraken2
        obs_reports, obs_outputs = _classify_kraken(manifest, common_args)

        self.assertIsInstance(obs_reports, Kraken2ReportDirectoryFormat)
        self.assertIsInstance(obs_outputs, Kraken2OutputDirectoryFormat)

        p1.assert_has_calls(
            [
                call(
                    cmd=[
                        "kraken2",
                        "--db",
                        "/some/where/db",
                        "--quick",
                        "--report",
                        exp_rep_fps[0],
                        "--output",
                        exp_out_fps[0],
                        "--use-names",
                        "samp1/bin1.fa",
                    ],
                    verbose=True,
                ),
                call(
                    cmd=[
                        "kraken2",
                        "--db",
                        "/some/where/db",
                        "--quick",
                        "--report",
                        exp_rep_fps[1],
                        "--output",
                        exp_out_fps[1],
                        "--use-names",
                        "samp1/bin2.fa",
                    ],
                    verbose=True,
                ),
                call(
                    cmd=[
                        "kraken2",
                        "--db",
                        "/some/where/db",
                        "--quick",
                        "--report",
                        exp_rep_fps[2],
                        "--output",
                        exp_out_fps[2],
                        "--use-names",
                        "samp2/bin1.fa",
                    ],
                    verbose=True,
                ),
            ]
        )
        p2.assert_has_calls(
            [
                call("samp1", "bin1", fake_output_dir, fake_report_dir),
                call("samp1", "bin2", fake_output_dir, fake_report_dir),
                call("samp2", "bin1", fake_output_dir, fake_report_dir),
            ]
        )
        p3.assert_has_calls(
            [
                call(("samp1", "bin1"), ANY, manifest.columns),
                call(("samp1", "bin2"), ANY, manifest.columns),
                call(("samp2", "bin1"), ANY, manifest.columns),
            ]
        )


if __name__ == "__main__":
    unittest.main()
