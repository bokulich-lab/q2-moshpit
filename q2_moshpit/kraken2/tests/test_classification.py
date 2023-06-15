# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import unittest
from subprocess import CalledProcessError

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
)
from qiime2 import Artifact

from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat, Kraken2DBDirectoryFormat,
)

from q2_types_genomics.per_sample_data._format import MultiMAGSequencesDirFmt
from unittest.mock import patch, ANY, call

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import moshpit

from q2_moshpit.kraken2.classification import (
    _get_seq_paths, _construct_output_paths, _classify_kraken2
)


class TestKraken2Classification(TestPluginBase):
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

    @patch("q2_moshpit.classification.Kraken2OutputDirectoryFormat")
    @patch("q2_moshpit.classification.Kraken2ReportDirectoryFormat")
    @patch("q2_moshpit.classification._get_seq_paths")
    @patch("q2_moshpit.classification._construct_output_paths")
    @patch("q2_moshpit.classification.run_command")
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
        obs_reports, obs_outputs = _classify_kraken2(manifest, common_args)

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

    @patch("q2_moshpit.classification.Kraken2OutputDirectoryFormat")
    @patch("q2_moshpit.classification.Kraken2ReportDirectoryFormat")
    @patch(
        "q2_moshpit.classification._get_seq_paths",
        return_value=(1, 2, [3])
    )
    @patch(
        "q2_moshpit.classification._construct_output_paths",
        return_value=(1, 2)
    )
    @patch("q2_moshpit.classification.run_command")
    def test_classify_kraken_exception(self, p1, p2, p3, p4, p5):
        manifest = MultiMAGSequencesDirFmt(
            self.get_data_path("mags"), "r"
        ).manifest.view(pd.DataFrame)
        common_args = ["--db", "/some/where/db", "--quick"]

        # run kraken2
        p1.side_effect = CalledProcessError(returncode=123, cmd="abc")
        with self.assertRaisesRegex(
            Exception,
            r'error was encountered .* \(return code 123\)'
        ):
            _classify_kraken2(manifest, common_args)

    @patch("q2_moshpit.classification._classify_kraken2")
    def test_classify_kraken_action(self, p1):
        seqs = Artifact.import_data(
            'SampleData[MAGs]', self.get_data_path("mags")
        )
        db = Artifact.import_data('Kraken2DB', self.get_data_path("db"))
        p1.return_value = (
            Kraken2ReportDirectoryFormat(
                self.get_data_path("reports-mags"), "r"
            ),
            Kraken2OutputDirectoryFormat(
                self.get_data_path("outputs-mags"), "r"
            ),
        )

        moshpit.actions.classify_kraken2(
            seqs=seqs, kraken2_db=db, threads=3, confidence=0.9, quick=True
        )

        exp_args = [
            '--threads', '3', '--confidence', '0.9',
            '--minimum-base-quality', '0', '--minimum-hit-groups', '2',
            '--quick', '--db', str(db.view(Kraken2DBDirectoryFormat).path)
        ]
        p1.assert_called_with(ANY, exp_args)


if __name__ == "__main__":
    unittest.main()
