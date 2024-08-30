# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import os
import unittest
from subprocess import CalledProcessError
from unittest.mock import patch, Mock, ANY

import numpy as np
import pandas as pd
from pandas._testing import assert_frame_equal
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt
)

from q2_moshpit.kaiju.classification import (
    _construct_feature_table, _rename_taxon, _encode_unclassified_ids,
    _fix_id_types, _process_kaiju_reports, _classify_kaiju
)
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import moshpit


class TestKaijuClassification(TestPluginBase):
    package = 'q2_moshpit.kaiju.tests'

    def setUp(self):
        super().setUp()
        self.classify_kaiju = moshpit.pipelines.classify_kaiju
        with open(self.get_data_path('taxa-map.json')) as f:
            self.taxa_map = json.load(f)

    def test_construct_feature_table(self):
        obs_table, obs_taxonomy = _construct_feature_table(
            table_fp=self.get_data_path('kaiju-table.tsv')
        )

        exp_table = pd.read_csv(
            self.get_data_path('kaiju-ft-ok.csv'), index_col=0
        )
        exp_table.columns.name = "taxon_id"
        exp_taxonomy = pd.read_csv(
            self.get_data_path('kaiju-taxonomy-ok.csv'), index_col=0
        )

        assert_frame_equal(exp_table, obs_table)
        assert_frame_equal(exp_taxonomy, obs_taxonomy)

    def test_rename_taxon_unspecified(self):
        obs_taxon = _rename_taxon("1236", self.taxa_map)
        exp_taxon = (
            "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
            "o__Unspecified;f__Unspecified;g__Unspecified;s__Unspecified"
        )

        self.assertEqual(exp_taxon, obs_taxon)

    def test_rename_taxon_full(self):
        obs_taxon = _rename_taxon("1773", self.taxa_map)
        exp_taxon = ("d__Bacteria;p__Actinobacteria;c__Actinomycetia;"
                     "o__Corynebacteriales;f__Mycobacteriaceae;"
                     "g__Mycobacterium;s__Mycobacterium tuberculosis")

        self.assertEqual(exp_taxon, obs_taxon)

    def test_rename_taxon_unclassified(self):
        obs_taxon = _rename_taxon("dW5jbGFz", self.taxa_map)
        exp_taxon = "d__unclassified"

        self.assertEqual(exp_taxon, obs_taxon)

    def test_encode_unclassified(self):
        df_in = pd.DataFrame({
            "taxon_name": ["Alpha", "Beta", "Alpha", "Gamma"],
            "taxon_id": [1, 2, 3, 4]
        })

        obs_table = _encode_unclassified_ids(df_in, "Alph")
        exp_table = pd.DataFrame({
            "taxon_name": ["Alpha", "Beta", "Alpha", "Gamma"],
            "taxon_id": ["QWxwaGE=", 2, "QWxwaGE=", 4]
        })

        assert_frame_equal(exp_table, obs_table)

    def test_fix_id_types(self):
        df_in = pd.DataFrame({
            "taxon_name": [
                "cannot be assigned", "Beta", "unclassified", "unclassified"
            ],
            "taxon_id": [np.nan, 21.0, np.nan, np.nan]
        })

        obs_table = _fix_id_types(df_in)
        exp_table = pd.DataFrame({
            "taxon_name": [
                "cannot be assigned", "Beta", "unclassified", "unclassified"
            ],
            "taxon_id": [
                "Y2Fubm90", 21, "dW5jbGFz", "dW5jbGFz"
            ]
        })

        assert_frame_equal(exp_table, obs_table, check_dtype=False)

    @patch("subprocess.run")
    @patch("q2_moshpit.kaiju.classification._construct_feature_table")
    def test_process_kaiju_reports_c_float(self, p1, p2):
        open(
            os.path.join(str(self.temp_dir.name), "sample1.out"), "w"
        ).close()
        open(
            os.path.join(str(self.temp_dir.name), "sample2.out"), "w"
        ).close()
        args = {
            "r": "species",
            "db": Mock(path=self.temp_dir.name),
            "exp": True,
            "u": True,
            "c": 0.6
        }

        _process_kaiju_reports(self.temp_dir.name, args)

        exp_cmd = [
            "kaiju2table", "-v", "-o",
            f"{self.temp_dir.name}/results.tsv",
            "-r", "species",
            "-t", f"{self.temp_dir.name}/nodes.dmp",
            "-n", f"{self.temp_dir.name}/names.dmp",
            "-l", "superkingdom,phylum,class,order,family,genus,species",
            "-e", "-u", "-m", "0.6",
            f"{self.temp_dir.name}/sample1.out",
            f"{self.temp_dir.name}/sample2.out"
        ]

        p1.assert_called_once_with(f"{self.temp_dir.name}/results.tsv")
        p2.assert_called_once_with(exp_cmd, check=True)

    @patch("subprocess.run")
    @patch("q2_moshpit.kaiju.classification._construct_feature_table")
    def test_process_kaiju_reports_c_int(self, p1, p2):
        open(
            os.path.join(str(self.temp_dir.name), "sample1.out"), "w"
        ).close()
        open(
            os.path.join(str(self.temp_dir.name), "sample2.out"), "w"
        ).close()
        args = {
            "r": "species",
            "db": Mock(path=self.temp_dir.name),
            "exp": True,
            "u": True,
            "c": 2
        }

        _process_kaiju_reports(self.temp_dir.name, args)

        exp_cmd = [
            "kaiju2table", "-v", "-o",
            f"{self.temp_dir.name}/results.tsv",
            "-r", "species",
            "-t", f"{self.temp_dir.name}/nodes.dmp",
            "-n", f"{self.temp_dir.name}/names.dmp",
            "-l", "superkingdom,phylum,class,order,family,genus,species",
            "-e", "-u", "-c", "2",
            f"{self.temp_dir.name}/sample1.out",
            f"{self.temp_dir.name}/sample2.out"
        ]

        p1.assert_called_once_with(f"{self.temp_dir.name}/results.tsv")
        p2.assert_called_once_with(exp_cmd, check=True)

    @patch("subprocess.run")
    @patch("q2_moshpit.kaiju.classification._process_kaiju_reports")
    def test_classify_kaiju_single(self, p1, p2):
        seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path("single-end"), "r"
        )
        db_path = self.temp_dir.name
        open(os.path.join(db_path, "kaiju_123.fmi"), "w").close()
        p1.return_value = [pd.DataFrame(), pd.DataFrame()]

        with patch("tempfile.TemporaryDirectory"):
            _classify_kaiju(
                seqs=seqs, db=Mock(path=self.temp_dir.name),
                z=3, a="greedy", e=2, m=10, s=66, evalue=0, x=False,
                r="class", c=0.1, exp=True, u=True
            )

        exp_cmd = [
            "kaiju-multi", "-v", "-z", "3", "-a", "greedy", "-e", "2",
            "-m", "10", "-s", "66", "-E", "0", "-X",
            "-t", os.path.join(db_path, "nodes.dmp"),
            "-f", os.path.join(db_path, "kaiju_123.fmi"),
            "-i", ANY, "-o", ANY
        ]
        obs_cmd = p2.call_args.args[0]

        p2.assert_called_once_with(exp_cmd, check=True)
        self.assertRegex(
            obs_cmd[-1], r".*sample1.out,.*sample2.out"
        )
        self.assertRegex(
            obs_cmd[-3], r".*reads1_R1.fastq.gz,.*reads2_R1.fastq.gz"
        )

    @patch("subprocess.run")
    @patch("q2_moshpit.kaiju.classification._process_kaiju_reports")
    def test_classify_kaiju_paired(self, p1, p2):
        seqs = SingleLanePerSamplePairedEndFastqDirFmt(
            self.get_data_path("paired-end"), "r"
        )
        db_path = self.temp_dir.name
        open(os.path.join(db_path, "kaiju_123.fmi"), "w").close()
        p1.return_value = [pd.DataFrame(), pd.DataFrame()]

        with patch("tempfile.TemporaryDirectory"):
            _classify_kaiju(
                seqs=seqs, db=Mock(path=self.temp_dir.name),
                z=3, a="greedy", e=2, m=10, s=66, evalue=0, x=True,
                r="class", c=0.1, exp=False, u=False
            )

        exp_cmd = [
            "kaiju-multi", "-v", "-z", "3", "-a", "greedy", "-e", "2",
            "-m", "10", "-s", "66", "-E", "0", "-x",
            "-t", os.path.join(db_path, "nodes.dmp"),
            "-f", os.path.join(db_path, "kaiju_123.fmi"),
            "-i", ANY, "-j", ANY, "-o", ANY
        ]
        obs_cmd = p2.call_args.args[0]

        p2.assert_called_once_with(exp_cmd, check=True)
        self.assertRegex(
            obs_cmd[-1], r".*sample1.out,.*sample2.out"
        )
        self.assertRegex(
            obs_cmd[-3], r".*reads1_R2.fastq.gz,.*reads2_R2.fastq.gz"
        )
        self.assertRegex(
            obs_cmd[-5], r".*reads1_R1.fastq.gz,.*reads2_R1.fastq.gz"
        )

    @patch("subprocess.run", side_effect=CalledProcessError(1, "hello"))
    def test_classify_kaiju_exception(self, p1):
        seqs = SingleLanePerSampleSingleEndFastqDirFmt(
            self.get_data_path("single-end"), "r"
        )
        open(os.path.join(self.temp_dir.name, "kaiju_123.fmi"), "w").close()

        with self.assertRaisesRegex(
                Exception, r"\(return code 1\), please inspect"
        ):
            _classify_kaiju(
                seqs=seqs, db=Mock(path=self.temp_dir.name),
                z=3, a="greedy", e=2, m=10, s=66, evalue=0, x=False,
                r="class", c=0.1, exp=True, u=True
            )


if __name__ == "__main__":
    unittest.main()
