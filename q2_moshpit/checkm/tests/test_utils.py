# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest

from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.checkm.utils import _get_plots_per_sample


class TestCheckMUtils(TestPluginBase):
    package = 'q2_moshpit.checkm.tests'

    def test_get_plots_per_sample(self):
        obs = _get_plots_per_sample({
            'plots_gc': {'samp1': 'abc', 'samp2': 'def'},
            'plots_nx': {'samp1': 'cba', 'samp2': 'fed'}
        })
        exp = {
            'samp1': {'plots_gc': 'abc', 'plots_nx': 'cba'},
            'samp2': {'plots_gc': 'def', 'plots_nx': 'fed'}
        }
        self.assertDictEqual(exp, obs)

    def test_get_plots_per_sample_uneven(self):
        with self.assertRaisesRegex(
            ValueError,
            r'.*Sample counts were: \[2, 1, 3\].'
        ):
            _get_plots_per_sample({
                'plots_gc': {'samp1': 'abc', 'samp2': 'def'},
                'plots_nx': {'samp1': 'cba'},
                'plots_coding': {'samp1': 'a', 'samp2': 'd', 'samp3': 'e'}
            })


if __name__ == '__main__':
    unittest.main()
