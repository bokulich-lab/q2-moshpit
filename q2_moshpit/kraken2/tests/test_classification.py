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

from q2_types_genomics.feature_data import MAGSequencesDirFmt
from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
    Kraken2DBDirectoryFormat,
)

from unittest.mock import patch, ANY, call

from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import moshpit

from q2_moshpit.kraken2.classification import (
    _get_seq_paths,
    _construct_output_paths,
    _classify_kraken2,
)


class TestKraken2Classification(TestPluginBase):
    package = "q2_moshpit.kraken2.tests"

    def setUp(self):
        super().setUp()

    def test_get_seq_paths_reads_single(self):
        manifest = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path("single-end"), mode="r"
        ).manifest.view(pd.DataFrame)

        expected = [(f"sample{i}", f"reads{i}_R1.fastq.gz") for i in range(1, 3)]
        for (index, row), (exp_sample, exp_fp) in zip(manifest.iterrows(), expected):
            _sample, fn = _get_seq_paths(index, row, manifest.columns)
            self.assertEqual(_sample, exp_sample)
            self.assertTrue(fn[0].endswith(exp_fp))

    def test_get_seq_paths_reads_paired(self):
        manifest = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path("paired-end"), mode="r"
        ).manifest.view(pd.DataFrame)

        expected = [
            (
                f"sample{i}",
                [f"reads{i}_R1.fastq.gz", f"reads{i}_R2.fastq.gz"],
            )
            for i in range(1, 3)
        ]
        for (index, row), (exp_sample, exp_fp) in zip(manifest.iterrows(), expected):
            _sample, fn = _get_seq_paths(index, row, manifest.columns)
            self.assertEqual(_sample, exp_sample)
            self.assertTrue(fn[0].endswith(exp_fp[0]))
            self.assertTrue(fn[1].endswith(exp_fp[1]))

    def test_construct_output_paths(self):
        _sample = "sample1"
        reports_dir = Kraken2ReportDirectoryFormat()
        outputs_dir = Kraken2OutputDirectoryFormat()

        exp_rep_fp = os.path.join(reports_dir.path, f"{_sample}.report.txt")
        exp_out_fp = os.path.join(outputs_dir.path, f"{_sample}.output.txt")
        obs_out_fp, obs_rep_fp = _construct_output_paths(
            _sample, outputs_dir, reports_dir
        )
        self.assertEqual(obs_rep_fp, exp_rep_fp)
        self.assertEqual(obs_out_fp, exp_out_fp)

    @patch("q2_moshpit.kraken2.classification.Kraken2OutputDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification.Kraken2ReportDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification._get_seq_paths")
    @patch("q2_moshpit.kraken2.classification._construct_output_paths")
    @patch("q2_moshpit.kraken2.classification.run_command")
    def test_classify_kraken_mags(self, p1, p2, p3, p4, p5):
        seqs = MAGSequencesDirFmt(self.get_data_path("mags-derep"), "r")
        common_args = ["--db", "/some/where/db", "--quick"]

        fake_report_dir = Kraken2ReportDirectoryFormat()
        fake_output_dir = Kraken2OutputDirectoryFormat()
        exp_out_fps = [
            os.path.join(
                fake_output_dir.path, "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.output.txt"
            ),
            os.path.join(
                fake_output_dir.path, "8894435a-c836-4c18-b475-8b38a9ab6c6b.output.txt"
            ),
        ]
        exp_rep_fps = [
            os.path.join(
                fake_report_dir.path, "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.report.txt"
            ),
            os.path.join(
                fake_report_dir.path, "8894435a-c836-4c18-b475-8b38a9ab6c6b.report.txt"
            ),
        ]

        p2.side_effect = list(zip(exp_out_fps, exp_rep_fps))
        p4.return_value = fake_report_dir
        p5.return_value = fake_output_dir

        # run kraken2
        obs_reports, obs_outputs = _classify_kraken2(seqs, common_args)

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
                        os.path.join(
                            seqs.path, "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.fasta"
                        ),
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
                        os.path.join(
                            seqs.path, "8894435a-c836-4c18-b475-8b38a9ab6c6b.fasta"
                        ),
                    ],
                    verbose=True,
                ),
            ]
        )
        p2.assert_has_calls(
            [
                call(
                    "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa",
                    fake_output_dir,
                    fake_report_dir,
                ),
                call(
                    "8894435a-c836-4c18-b475-8b38a9ab6c6b",
                    fake_output_dir,
                    fake_report_dir,
                ),
            ]
        )
        p3.assert_not_called()

    @patch("q2_moshpit.kraken2.classification.Kraken2OutputDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification.Kraken2ReportDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification._get_seq_paths")
    @patch("q2_moshpit.kraken2.classification._construct_output_paths")
    @patch("q2_moshpit.kraken2.classification.run_command")
    def test_classify_kraken_reads(self, p1, p2, p3, p4, p5):
        seqs = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path("paired-end"), "r"
        )
        manifest = seqs.manifest.view(pd.DataFrame)
        common_args = ["--db", "/some/where/db", "--quick"]

        fake_report_dir = Kraken2ReportDirectoryFormat()
        fake_output_dir = Kraken2OutputDirectoryFormat()
        exp_out_fps = [
            os.path.join(fake_output_dir.path, "sample1.output.txt"),
            os.path.join(fake_output_dir.path, "sample2.output.txt"),
        ]
        exp_rep_fps = [
            os.path.join(fake_report_dir.path, "sample1.report.txt"),
            os.path.join(fake_report_dir.path, "sample2.report.txt"),
        ]

        p2.side_effect = list(zip(exp_out_fps, exp_rep_fps))
        p3.side_effect = [
            ("sample1", ["reads1_R1.fastq.gz", "reads1_R2.fastq.gz"]),
            ("sample2", ["reads2_R1.fastq.gz", "reads2_R2.fastq.gz"]),
        ]
        p4.return_value = fake_report_dir
        p5.return_value = fake_output_dir

        # run kraken2
        obs_reports, obs_outputs = _classify_kraken2(seqs, common_args)

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
                        "--paired",
                        "--report",
                        exp_rep_fps[0],
                        "--output",
                        exp_out_fps[0],
                        "reads1_R1.fastq.gz",
                        "reads1_R2.fastq.gz",
                    ],
                    verbose=True,
                ),
                call(
                    cmd=[
                        "kraken2",
                        "--db",
                        "/some/where/db",
                        "--quick",
                        "--paired",
                        "--report",
                        exp_rep_fps[1],
                        "--output",
                        exp_out_fps[1],
                        "reads2_R1.fastq.gz",
                        "reads2_R2.fastq.gz",
                    ],
                    verbose=True,
                ),
            ]
        )
        p2.assert_has_calls(
            [
                call("sample1", fake_output_dir, fake_report_dir),
                call("sample2", fake_output_dir, fake_report_dir),
            ]
        )
        p3.assert_has_calls(
            [
                call("sample1", ANY, list(manifest.columns)),
                call("sample2", ANY, list(manifest.columns)),
            ]
        )

    @patch("q2_moshpit.kraken2.classification.Kraken2OutputDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification.Kraken2ReportDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification._get_seq_paths", return_value=(1, 2, [3]))
    @patch(
        "q2_moshpit.kraken2.classification._construct_output_paths", return_value=(1, 2)
    )
    @patch("q2_moshpit.kraken2.classification.run_command")
    def test_classify_kraken_exception(self, p1, p2, p3, p4, p5):
        seqs = MAGSequencesDirFmt(self.get_data_path("mags-derep"), "r")
        common_args = ["--db", "/some/where/db", "--quick"]

        # run kraken2
        p1.side_effect = CalledProcessError(returncode=123, cmd="abc")
        with self.assertRaisesRegex(
            Exception, r"error was encountered .* \(return code 123\)"
        ):
            _classify_kraken2(seqs, common_args)

    @patch("q2_moshpit.kraken2.classification._classify_kraken2")
    def test_classify_kraken_action(self, p1):
        seqs = Artifact.import_data(
            "FeatureData[MAG]", self.get_data_path("mags-derep")
        )
        db = Artifact.import_data("Kraken2DB", self.get_data_path("db"))
        p1.return_value = (
            Kraken2ReportDirectoryFormat(self.get_data_path("reports-mags"), "r"),
            Kraken2OutputDirectoryFormat(self.get_data_path("outputs-mags"), "r"),
        )

        moshpit.actions.classify_kraken2(
            seqs=seqs, kraken2_db=db, threads=3, confidence=0.9, quick=True
        )

        exp_args = [
            "--threads",
            "3",
            "--confidence",
            "0.9",
            "--minimum-base-quality",
            "0",
            "--minimum-hit-groups",
            "2",
            "--quick",
            "--db",
            str(db.view(Kraken2DBDirectoryFormat).path),
        ]
        p1.assert_called_with(ANY, exp_args)


if __name__ == "__main__":
    unittest.main()
