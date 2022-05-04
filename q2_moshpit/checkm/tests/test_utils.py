# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import pandas as pd
from pandas._testing import assert_frame_equal
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.checkm.checkm import (
    _parse_single_checkm_report, _parse_checkm_reports
)
from q2_moshpit.checkm.utils import (
    _extract_checkm_stats, BinStatistics, MarkerCounts, SampleStatistics,
    _extract_pattern
)


class TestCheckMUtils(TestPluginBase):
    package = 'q2_moshpit.checkm.tests'

    def setUp(self):
        super().setUp()
        self.report_line1 = \
            '  bin1      g__Mycobacterium (UID1816)       100         693   ' \
            '         300        0    693   0   0   0   0       100.00      ' \
            '     0.00               0.00          '
        self.report_line2 = \
            '  bin_1      g__Strain_ales (UID18)       0         693   '\
            '         30        10    682   1   0   1   0       0.14      ' \
            '     99.19               15.84'
        self.sample_stats = SampleStatistics(
            sample_id='samp1',
            bins=[
                BinStatistics(
                    bin_id='bin1',
                    marker_lineage='g__Mycobacterium',
                    marker_lineage_uid='UID1816',
                    genomes=100, markers=693, marker_sets=300,
                    marker_counts=MarkerCounts(0, 693, 0, 0, 0, 0),
                    completeness=100.0, contamination=0.0,
                    strain_heterogeneity=0.0
                ),
                BinStatistics(
                    bin_id='bin2',
                    marker_lineage='o__Pseudomonadales',
                    marker_lineage_uid='UID4488',
                    genomes=185, markers=813, marker_sets=308,
                    marker_counts=MarkerCounts(97, 707, 9, 0, 0, 0),
                    completeness=93.87, contamination=1.31,
                    strain_heterogeneity=0.0
                )
            ]
        )

    def test_extract_pattern_ints(self):
        obs = _extract_pattern(r'\s+(\d{1,})', self.report_line1, int)
        exp = [100, 693, 300, 0, 693, 0, 0, 0, 0, 100, 0, 0]
        self.assertListEqual(exp, obs)

    def test_extract_pattern_ints_les_zeros(self):
        obs = _extract_pattern(r'\s+(\d{1,})', self.report_line2, int)
        exp = [0, 693, 30, 10, 682, 1, 0, 1, 0, 0, 99, 15]
        self.assertListEqual(exp, obs)

    def test_extract_pattern_floats(self):
        obs = _extract_pattern(r'\d{1,3}\.\d{2}\s*', self.report_line1, float)
        exp = [100.0, 0.0, 0.0]
        self.assertListEqual(exp, obs)

    def test_extract_pattern_floats_less_zeros(self):
        obs = _extract_pattern(r'\d{1,3}\.\d{2}\s*', self.report_line2, float)
        exp = [0.14, 99.19, 15.84]
        self.assertListEqual(exp, obs)

    def test_extract_pattern_words(self):
        obs = _extract_pattern(
            r'(\S+)\s+(\D__\w+)\s\((UID\d+)\)', self.report_line1
        )
        exp = [('bin1', 'g__Mycobacterium', 'UID1816')]
        self.assertListEqual(exp, obs)

    def test_extract_pattern_words_more_underscores(self):
        obs = _extract_pattern(
            r'(\S+)\s+(\D__\w+)\s\((UID\d+)\)', self.report_line2
        )
        exp = [('bin_1', 'g__Strain_ales', 'UID18')]
        self.assertListEqual(exp, obs)

    def test_extract_all(self):
        obs = _extract_checkm_stats(self.report_line1)
        exp = BinStatistics(
            bin_id='bin1', marker_lineage='g__Mycobacterium',
            marker_lineage_uid='UID1816',
            genomes=100, markers=693, marker_sets=300,
            marker_counts=MarkerCounts(0, 693, 0, 0, 0, 0),
            completeness=100.0, contamination=0.0, strain_heterogeneity=0.0
        )
        self.assertEqual(exp, obs)

    def test_sample_stats_to_df(self):
        obs = self.sample_stats.to_df()
        exp = pd.read_csv(
            self.get_data_path('checkm_report_df1.tsv'),
            sep='\t', index_col=None
        )
        assert_frame_equal(exp, obs)

    def test_parse_single_checkm_report(self):
        obs = _parse_single_checkm_report(
            'samp1', self.get_data_path('raw_report1.txt')
        )
        exp = pd.read_csv(
            self.get_data_path('checkm_report_df1.tsv'),
            sep='\t', index_col=None
        )
        assert_frame_equal(exp, obs)

    def test_parse_single_checkm_report_with_error(self):
        with self.assertRaisesRegexp(
            ValueError,
            r'counts \(694\) is different from .* marker count \(693\)'
        ):
            _parse_single_checkm_report(
                'samp1', self.get_data_path('raw_report1_wrong.txt')
            )

    def test_parse_multiple_checkm_report(self):
        obs = _parse_checkm_reports(
            {
                'samp1': self.get_data_path('raw_report1.txt'),
                'samp2': self.get_data_path('raw_report2.txt')
            }
        )
        exp = pd.read_csv(
            self.get_data_path('checkm_report_df_all.tsv'),
            sep='\t', index_col=None
        )
        assert_frame_equal(exp, obs)


if __name__ == '__main__':
    unittest.main()
