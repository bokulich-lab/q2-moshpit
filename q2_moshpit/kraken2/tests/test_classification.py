# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from pathlib import Path
import unittest
from subprocess import CalledProcessError

import pandas as pd
from unittest.mock import patch, ANY, call

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
)
from q2_types_genomics.feature_data import MAGSequencesDirFmt
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat, Kraken2DBDirectoryFormat,
)
from q2_moshpit.kraken2.classification import (
    _get_seq_paths, _construct_output_paths, _classify_kraken2,
    classify_kraken2
)

from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import moshpit


class TestClassifyKraken2Helpers(TestPluginBase):
    package = "q2_moshpit.kraken2.tests"

    def setUp(self):
        super().setUp()

    def test_get_seq_paths_reads_single(self):
        manifest = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path("single-end"), mode="r"
        ).manifest.view(pd.DataFrame)

        expected = [
            (f"sample{i}", f"reads{i}_R1.fastq.gz")
            for i in range(1, 3)
        ]
        for ((index, row), (exp_sample, exp_fp)) in zip(
            manifest.iterrows(), expected
        ):
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
        for ((index, row), (exp_sample, exp_fp)) in zip(
            manifest.iterrows(), expected
        ):
            _sample, fn = _get_seq_paths(index, row, manifest.columns)
            self.assertEqual(_sample, exp_sample)
            self.assertTrue(fn[0].endswith(exp_fp[0]))
            self.assertTrue(fn[1].endswith(exp_fp[1]))

    def test_construct_output_paths(self):
        _sample = "sample1"
        reports_dir = Kraken2ReportDirectoryFormat()
        outputs_dir = Kraken2OutputDirectoryFormat()

        exp_rep_fp = os.path.join(
            reports_dir.path, f"{_sample}.report.txt"
        )
        exp_out_fp = os.path.join(
            outputs_dir.path, f"{_sample}.output.txt"
        )
        obs_out_fp, obs_rep_fp = _construct_output_paths(
            _sample, outputs_dir, reports_dir
        )
        self.assertEqual(obs_rep_fp, exp_rep_fp)
        self.assertEqual(obs_out_fp, exp_out_fp)


class TestClassifyKraken2HasCorrectCalls(TestPluginBase):
    package = "q2_moshpit.kraken2.tests"

    def setUp(self):
        super().setUp()

    @patch("q2_moshpit.kraken2.classification.Kraken2OutputDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification.Kraken2ReportDirectoryFormat")
    @patch(
        "q2_moshpit.kraken2.classification._get_seq_paths",
        return_value=(1, 2, [3])
    )
    @patch(
        "q2_moshpit.kraken2.classification._construct_output_paths",
        return_value=(1, 2)
    )
    @patch("q2_moshpit.kraken2.classification.run_command")
    def test_exception(self, p1, p2, p3, p4, p5):
        seqs = MAGSequencesDirFmt(self.get_data_path("mags-derep"), "r")
        common_args = ["--db", "/some/where/db", "--quick"]

        # run kraken2
        p1.side_effect = CalledProcessError(returncode=123, cmd="abc")
        with self.assertRaisesRegex(
            Exception,
            r'error was encountered .* \(return code 123\)'
        ):
            _classify_kraken2(seqs, common_args)

    @patch("q2_moshpit.kraken2.classification._classify_kraken2")
    def test_action(self, p1):
        seqs = Artifact.import_data(
            'FeatureData[MAG]', self.get_data_path("mags-derep")
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

    @patch("q2_moshpit.kraken2.classification.Kraken2OutputDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification.Kraken2ReportDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification._get_seq_paths")
    @patch("q2_moshpit.kraken2.classification._construct_output_paths")
    @patch("q2_moshpit.kraken2.classification.run_command")
    def test_reads(
        self, p1, p2, p3, p4, p5
    ):
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
    @patch("q2_moshpit.kraken2.classification._get_seq_paths")
    @patch("q2_moshpit.kraken2.classification.run_command")
    def test_contigs(
        self,
        run_command_mock,
        _get_seq_paths_mock,
        report_format_mock,
        output_format_mock
    ):
        samples_dir = self.get_data_path(os.path.join('simulated-sequences', 'contigs'))
        contigs = ContigSequencesDirFmt(samples_dir, "r")

        common_args = ["--db", "/some/where/db", "--quick"]

        fake_output_dir = Kraken2OutputDirectoryFormat()
        fake_report_dir = Kraken2ReportDirectoryFormat()

        samples = ('ba', 'mm', 'sa', 'se')
        exp_output_fps = []
        exp_report_fps = []
        for sample in samples:
            exp_output_fps.append(
                os.path.join(fake_output_dir.path, f'{sample}.output.txt')
            )
            exp_report_fps.append(
                os.path.join(fake_report_dir.path, f'{sample}.report.txt')
            )

        output_format_mock.return_value = fake_output_dir
        report_format_mock.return_value = fake_report_dir

        obs_reports, obs_outputs = _classify_kraken2(contigs, common_args)
        self.assertIsInstance(obs_reports, Kraken2ReportDirectoryFormat)
        self.assertIsInstance(obs_outputs, Kraken2OutputDirectoryFormat)

        calls = []
        for i, sample in enumerate(samples):
            calls.append(call(
                cmd=[
                    "kraken2",
                    "--db",
                    "/some/where/db",
                    "--quick",
                    "--report",
                    exp_report_fps[i],
                    "--output",
                    exp_output_fps[i],
                    os.path.join(
                        contigs.path,
                        f'{sample}_contigs.fasta'
                    )
                ],
                verbose=True
            ))
        run_command_mock.assert_has_calls(calls, any_order=True)

        _get_seq_paths_mock.assert_not_called()

    @patch("q2_moshpit.kraken2.classification.Kraken2OutputDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification.Kraken2ReportDirectoryFormat")
    @patch("q2_moshpit.kraken2.classification._get_seq_paths")
    @patch("q2_moshpit.kraken2.classification._construct_output_paths")
    @patch("q2_moshpit.kraken2.classification.run_command")
    def test_mags(self, p1, p2, p3, p4, p5):
        seqs = MAGSequencesDirFmt(self.get_data_path("mags-derep"), "r")
        common_args = ["--db", "/some/where/db", "--quick"]

        fake_report_dir = Kraken2ReportDirectoryFormat()
        fake_output_dir = Kraken2OutputDirectoryFormat()
        exp_out_fps = [
            os.path.join(
                fake_output_dir.path,
                "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.output.txt"
            ),
            os.path.join(
                fake_output_dir.path,
                "8894435a-c836-4c18-b475-8b38a9ab6c6b.output.txt"
            ),
        ]
        exp_rep_fps = [
            os.path.join(
                fake_report_dir.path,
                "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.report.txt"
            ),
            os.path.join(
                fake_report_dir.path,
                "8894435a-c836-4c18-b475-8b38a9ab6c6b.report.txt"
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
                            seqs.path,
                            "3b72d1a7-ddb0-4dc7-ac36-080ceda04aaa.fasta"
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
                            seqs.path,
                            "8894435a-c836-4c18-b475-8b38a9ab6c6b.fasta"
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
                    fake_output_dir, fake_report_dir
                ),
                call(
                    "8894435a-c836-4c18-b475-8b38a9ab6c6b",
                    fake_output_dir, fake_report_dir
                ),
            ]
        )
        p3.assert_not_called()


class TestClassifyKraken2Reads(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        datadir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'data'
        )

        db_path = os.path.join(datadir, 'simulated-sequences', 'kraken2-db')
        reads_path = os.path.join(datadir, 'simulated-sequences', 'reads')

        db = Kraken2DBDirectoryFormat(db_path, 'r')
        samples = SingleLanePerSamplePairedEndFastqDirFmt(reads_path, 'r')

        cls.reports, cls.outputs = classify_kraken2(samples, db)
        cls.output_views = cls.outputs.reports.iter_views(pd.DataFrame)
        cls.report_views = cls.reports.reports.iter_views(pd.DataFrame)

        cls.sample_id_to_ncbi_id = {
            'ba': {1392},   # bacillus anthracis
            'mm': {10090},  # mus musculus
            'sa': {1280},   # staph aureus
            'se': {1282},   # staph epidermidis
            'ba-mm-mixed': {1392, 10090}
        }

    def test_formats(self):
        self.assertIsInstance(self.reports, Kraken2ReportDirectoryFormat)
        self.assertIsInstance(self.outputs, Kraken2OutputDirectoryFormat)
        self.reports.validate()
        self.outputs.validate()

    def test_reads(self):
        samples_of_interest = ('ba', 'mm', 'sa', 'se', 'ba-mm-mixed')

        def filter_views(arg):
            path, _ = arg
            return Path(path.stem).stem in samples_of_interest

        output_views = filter(filter_views, self.output_views)
        report_views = filter(filter_views, self.report_views)

        for path, df in output_views:
            sample_id = str(path).rsplit('.output.txt')[0]

            # the expected number of records are in the output
            self.assertEqual(len(df), 25)

            # all reads are classified
            self.assertEqual({'C'}, set(df['classification']))

            # all reads are classified correctly
            self.assertEqual(
                set(df['taxon_id']),
                self.sample_id_to_ncbi_id[sample_id]
            )

        for path, df in report_views:
            sample_id = str(path).rsplit('.report.txt')[0]

            # the dataframe is non-empty
            self.assertGreater(len(df), 0)

            # the correct taxonomy id(s) is present somewhere in the
            # classification tree, and none of the others are present
            exp = self.sample_id_to_ncbi_id[sample_id]
            obs = set(df['taxon_id'])
            all_samples = set().union(
                *[s for _, s in self.sample_id_to_ncbi_id.items()]
            )
            exp_missing = all_samples - exp
            self.assertEqual(exp & obs, exp)
            self.assertFalse(exp_missing & obs)

    def test_nonsense_reads(self):
        samples_of_interest = ('nonsense')

        def filter_views(arg):
            path, _ = arg
            return Path(path.stem).stem in samples_of_interest

        output_views = filter(filter_views, self.output_views)
        report_views = filter(filter_views, self.report_views)

        _, df = list(output_views)[0]

        # the expected number of records are in the output
        self.assertEqual(len(df), 25)

        # the sequences are unclassified
        self.assertEqual({'U'}, set(df['classification']))

        _, df = list(report_views)[0]

        # the reports file has one line for all unclassified sequences
        self.assertEqual(len(df), 1)

        # none of the db taxonomy ids are present in the report
        exp = {0}
        obs = set(df['taxon_id'])
        self.assertEqual(exp, obs)

    # TODO: need to decide what to do here, currently empty report files
    # raise a pandas EmptyDataError that makes validation fail
    # also, kraken2 doesnt output the output.txt file for empty inputs...
    # probably need to just disallow any empty input files
    def test_empty_reads(self):
        pass


class TestClassifyKraken2Contigs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        datadir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'data'
        )

        db_path = os.path.join(datadir, 'simulated-sequences', 'kraken2-db')
        contigs_path = os.path.join(datadir, 'simulated-sequences', 'contigs')

        db = Kraken2DBDirectoryFormat(db_path, 'r')
        samples = ContigSequencesDirFmt(contigs_path, 'r')

        cls.reports, cls.outputs = classify_kraken2(samples, db)
        cls.output_views = cls.outputs.reports.iter_views(pd.DataFrame)
        cls.report_views = cls.reports.reports.iter_views(pd.DataFrame)

        cls.sample_id_to_ncbi_id = {
            'ba': {1392},   # bacillus anthracis
            'mm': {10090},  # mus musculus
            'sa': {1280},   # staph aureus
            'se': {1282},    # staph epidermidis
            'ba-mm-mixed': {1392, 10090}
        }

    def test_formats(self):
        self.assertIsInstance(self.reports, Kraken2ReportDirectoryFormat)
        self.assertIsInstance(self.outputs, Kraken2OutputDirectoryFormat)
        self.reports.validate()
        self.outputs.validate()

    def test_contigs(self):
        samples_of_interest = ('ba', 'mm', 'sa', 'se', 'ba-mm-mixed')

        def filter_views(arg):
            path, _ = arg
            return Path(path.stem).stem in samples_of_interest

        output_views = filter(filter_views, self.output_views)
        report_views = filter(filter_views, self.report_views)

        for path, df in output_views:
            sample_id = str(path).rsplit('.output.txt')[0]

            # the expected number of records are in the output
            self.assertEqual(len(df), 20)

            # all contigs are classified
            self.assertEqual({'C'}, set(df['classification']))

            # all contigs are classified correctly
            self.assertEqual(
                set(df['taxon_id']),
                self.sample_id_to_ncbi_id[sample_id]
            )

        for path, df in report_views:
            sample_id = str(path).rsplit('.report.txt')[0]

            # the dataframe is non-empty
            self.assertGreater(len(df), 0)

            # the correct taxonomy id(s) is present somewhere in the
            # classification tree, and none of the others are present
            exp = self.sample_id_to_ncbi_id[sample_id]
            obs = set(df['taxon_id'])
            all_samples = set().union(
                *[s for _, s in self.sample_id_to_ncbi_id.items()]
            )
            exp_missing = all_samples - exp
            self.assertEqual(exp & obs, exp)
            self.assertFalse(exp_missing & obs)


class TestClassifyKraken2MAGs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        datadir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), 'data'
        )

        db_path = os.path.join(datadir, 'simulated-sequences', 'kraken2-db')
        mags_path = os.path.join(datadir, 'simulated-sequences', 'mags')

        db = Kraken2DBDirectoryFormat(db_path, 'r')
        samples = MAGSequencesDirFmt(mags_path, 'r')

        cls.reports, cls.outputs = classify_kraken2(samples, db)
        cls.output_views = cls.outputs.reports.iter_views(pd.DataFrame)
        cls.report_views = cls.reports.reports.iter_views(pd.DataFrame)

        cls.uuid_to_sample = {
            '9231448e-b591-4afc-9d8a-5255b1a24f08': 'ba',
            '7797bbd1-4f3c-4482-9828-fa4be13c9977': 'mm',
            '5693d0e1-be8e-40ab-9427-94a0ffc62963': 'sa',
            '8adb2c2f-bb49-4b1a-a9ac-daf985b35070': 'se',
        }
        cls.sample_id_to_ncbi_id = {
            'ba': {1392},   # bacillus anthracis
            'mm': {10090},  # mus musculus
            'sa': {1280},   # staph aureus
            'se': {1282},   # staph epidermidis
        }

    def test_formats(self):
        self.assertIsInstance(self.reports, Kraken2ReportDirectoryFormat)
        self.assertIsInstance(self.outputs, Kraken2OutputDirectoryFormat)
        self.reports.validate()
        self.outputs.validate()

    def test_mags(self):
        for path, df in self.output_views:
            mag_id = str(path).rsplit('.output.txt')[0]
            sample_id = self.uuid_to_sample[mag_id]

            # the expected number of records are in the output
            self.assertGreater(len(df), 1)

            # all mags are classified
            self.assertEqual({'C'}, set(df['classification']))

            # all mags are classified correctly
            self.assertEqual(
                set(df['taxon_id']),
                self.sample_id_to_ncbi_id[sample_id]
            )

        for path, df in self.report_views:
            mag_id = str(path).rsplit('.report.txt')[0]
            sample_id = self.uuid_to_sample[mag_id]

            # the dataframe is non-empty
            self.assertGreater(len(df), 0)

            # the correct taxonomy id(s) is present somewhere in the
            # classification tree, and none of the others are present
            exp = self.sample_id_to_ncbi_id[sample_id]
            obs = set(df['taxon_id'])
            all_samples = set().union(
                *[s for _, s in self.sample_id_to_ncbi_id.items()]
            )
            exp_missing = all_samples - exp
            self.assertEqual(exp & obs, exp)
            self.assertFalse(exp_missing & obs)


if __name__ == "__main__":
    unittest.main()
