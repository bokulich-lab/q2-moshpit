# ----------------------------------------------------------------------------
# Copyright (c) 2023-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import shutil
import tempfile

import pandas as pd
from pandas._testing import assert_frame_equal
from q2_moshpit.kraken2 import kraken2_to_features, kraken2_to_mag_features
from qiime2.plugin.testing import TestPluginBase

from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat,
)


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestKrakenSelect(TestPluginBase):
    package = "q2_moshpit.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.temp_dir = tempfile.mkdtemp()

        fp = self.get_data_path(
            "kraken2-reports-select/kraken2_reads_table.json"
        )
        self.kraken2_reads_table = pd.read_json(fp)
        self.kraken2_reads_table.columns = [
            str(x) for x in self.kraken2_reads_table.columns
        ]

        fp = self.get_data_path(
            "kraken2-reports-select/kraken2_reads_table_filtered.json"
        )
        self.kraken2_reads_table_filtered = pd.read_json(fp)
        self.kraken2_reads_table_filtered.columns = [
            str(x) for x in self.kraken2_reads_table_filtered.columns
        ]

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_kraken2_to_features_coverage_default(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("kraken2-reports-select/samples"), "r"
        )
        obs_table, obs_taxonomy = kraken2_to_features(reports)

        assert_frame_equal(obs_table, self.kraken2_reads_table_filtered)
        # TODO: test taxonomy

    def test_kraken2_to_features_coverage_0(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("kraken2-reports-select/samples"), "r"
        )
        obs_table, obs_taxonomy = kraken2_to_features(reports, 0.0)

        assert_frame_equal(obs_table, self.kraken2_reads_table)
        # TODO: test taxonomy

    def test_kraken2_to_mag_features_default(self):
        reports = Kraken2ReportDirectoryFormat(
            self.get_data_path("reports-mags"), "r"
        )
        hits = Kraken2OutputDirectoryFormat(
            self.get_data_path("outputs-mags"), "r"
        )
        obs_table, obs_taxonomy = kraken2_to_mag_features(
            reports, hits, 0.0
        )
