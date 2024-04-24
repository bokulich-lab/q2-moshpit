# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.abundance.abundance import rpkm, tpm


class TestAbundance(TestPluginBase):
    package = "q2_moshpit.abundance.tests"

    def setUp(self):
        super().setUp()
        self.df = pd.DataFrame({
            'sample_id': [*['s1']*4, *['s2']*4, *['s3']*4],
            'mag_id': ['m1', 'm2', 'm3', 'm4']*3,
            'length': [2000, 4000, 1000, 10000]*3,
            'numreads': [10, 20, 5, 0, 12, 25, 8, 0, 30, 60, 15, 1],
        })

    def test_rpkm(self):
        obs = rpkm(self.df)
        exp = pd.Series(
            [142857.14, 142857.14, 142857.14, 0.0,
             133333.33, 138888.89, 177777.78, 0.0,
             141509.43, 141509.43, 141509.43, 943.40]
        )
        pd.testing.assert_series_equal(obs, exp)

    def test_tpm(self):
        obs = tpm(self.df)
        exp = pd.Series(
            [333333.33, 333333.33, 333333.33, 0.0,
             296296.3, 308641.98, 395061.73, 0.0,
             332594.24, 332594.24, 332594.24, 2217.29]
        )
        pd.testing.assert_series_equal(obs, exp)
