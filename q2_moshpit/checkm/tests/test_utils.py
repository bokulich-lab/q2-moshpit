# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
import unittest

import pandas as pd
from pandas._testing import assert_frame_equal
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.checkm.checkm import (
    _parse_single_checkm_report, _parse_checkm_reports
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

    def read_in_checkm_report(self, fp: str) -> pd.DataFrame:
        df = pd.read_csv(self.get_data_path(fp), sep='\t', index_col=None)
        for col in df.columns:
            if col.startswith('gcn'):
                df[col] = df[col].apply(json.loads)
        return df

    def test_parse_single_checkm_report(self):
        obs = _parse_single_checkm_report(
            'samp1', self.get_data_path('bin_stats_ext1.tsv')
        )
        exp = self.read_in_checkm_report('checkm_report_df1.tsv')

        assert_frame_equal(exp, obs, check_less_precise=2)

    def test_parse_multiple_checkm_reports(self):
        obs = _parse_checkm_reports(
            {
                'samp1': self.get_data_path('bin_stats_ext1.tsv'),
                'samp2': self.get_data_path('bin_stats_ext2.tsv')
            }
        )
        exp = self.read_in_checkm_report('checkm_report_df_all.tsv')

        assert_frame_equal(exp, obs, check_less_precise=2)


if __name__ == '__main__':
    unittest.main()
