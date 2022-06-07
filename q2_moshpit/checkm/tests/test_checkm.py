# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import contextlib
import json
import os
import shutil
import tempfile
import unittest
from unittest.mock import patch, call
from zipfile import ZipFile

import pandas as pd
from pandas._testing import assert_frame_equal
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.checkm.checkm import (
    _parse_single_checkm_report, _parse_checkm_reports, _classify_completeness,
    _draw_checkm_plots, _zip_checkm_plots, _evaluate_bins
)
from q2_moshpit.checkm.utils import _get_plots_per_sample
from q2_types_genomics.per_sample_data import MultiMAGSequencesDirFmt


class TestCheckM(TestPluginBase):
    package = 'q2_moshpit.checkm.tests'

    def setUp(self):
        super().setUp()
        with contextlib.ExitStack() as stack:
            self._tmp = stack.enter_context(tempfile.TemporaryDirectory())
            self.addCleanup(stack.pop_all().close)
        self.bins = MultiMAGSequencesDirFmt(self.get_data_path('bins'), 'r')
        self.db_path = 'some/where/else'

    def read_in_checkm_report(self, fp: str) -> pd.DataFrame:
        df = pd.read_csv(self.get_data_path(fp), sep='\t', index_col=None)
        for col in df.columns:
            if col.startswith('gcn'):
                df[col] = df[col].apply(json.loads)
        return df

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

    def test_classify_completeness(self):
        self.assertEqual('near', _classify_completeness(90.5))
        self.assertEqual('substantial', _classify_completeness(75.0))
        self.assertEqual('moderate', _classify_completeness(52.0))
        self.assertEqual('partial', _classify_completeness(25.0))

    @patch('subprocess.run')
    def test_draw_checkm_plots(self, p1):
        obs_plots = _draw_checkm_plots(
            results_dir=self._tmp, bins=self.bins,
            db_path=self.db_path, plot_type='gc'
        )

        exp_cmds = [[
            'checkm', 'gc_plot', '-x', 'fasta', '--image_type', 'svg',
            '--font_size', '10', self.get_data_path(f'bins/samp{x}'),
            os.path.join(self._tmp, 'plots', 'gc', f'samp{x}'),
            '50', '75', '90'
        ] for x in range(1, 3)]
        exp_calls = [
            call(cmd, env={**os.environ, 'CHECKM_DATA_PATH': self.db_path},
                 check=True)
            for cmd in exp_cmds
        ]
        exp_plots = {
            f'samp{x}': os.path.join(self._tmp, 'plots', 'gc', f'samp{x}')
            for x in range(1, 3)
        }

        p1.assert_has_calls(exp_calls)
        self.assertDictEqual(obs_plots, exp_plots)

    @patch.object(ZipFile, 'write')
    def test_zip_checkm_plots(self, p1):
        fake_archive = os.path.join(self._tmp, 'plots.zip')
        fake_plots = {
            'samp1': {
                'gc': self.get_data_path('plots/gc/samp1'),
                'nx': self.get_data_path('plots/nx/samp1')
            },
            'samp2': {
                'gc': self.get_data_path('plots/gc/samp2'),
                'nx': self.get_data_path('plots/nx/samp2')
            }
        }

        _zip_checkm_plots(fake_plots, fake_archive)

        exp_calls = [
            call(
                os.path.join(self.get_data_path('plots'), x, y),
                arcname=os.path.join(x, y)
            ) for x, y in [
                ('gc/samp1/', 'gc.plot1.svg'), ('nx/samp1/', 'nx.plot1.svg'),
                ('gc/samp2/', 'gc.plot1.svg'), ('gc/samp2/', 'gc.plot2.svg'),
                ('nx/samp2/', 'nx.plot1.svg'), ('nx/samp2/', 'nx.plot2.svg')
            ]
        ]
        p1.assert_has_calls(exp_calls, any_order=True)

    @patch('subprocess.run')
    def test_evaluate_bins(self, p1):
        shutil.copytree(
            self.get_data_path('checkm_reports'), self._tmp, dirs_exist_ok=True
        )
        obs_fps = _evaluate_bins(
            results_dir=self._tmp, bins=self.bins, db_path=self.db_path,
            common_args=['--reduced_tree', '--threads', '2']
        )

        base_cmd = [
            'checkm', 'lineage_wf', '--reduced_tree', '--threads', '2',
            '-x', 'fasta'
        ]
        exp_calls = [
            call(
                [*base_cmd, self.get_data_path(f'bins/samp{x}'),
                 os.path.join(self._tmp,  f'samp{x}')],
                env={**os.environ, 'CHECKM_DATA_PATH': self.db_path},
                check=True
            ) for x in range(1, 3)
        ]
        p1.assert_has_calls(exp_calls)
        self.assertDictEqual(
            obs_fps,
            {f'samp{x}': os.path.join(
                self._tmp, f'samp{x}', 'storage', 'bin_stats_ext.tsv'
            ) for x in range(1, 3)}
        )

    @patch('subprocess.run')
    def test_evaluate_bins_report_missing(self, p1):
        missing_report_fp = os.path.join(
            self._tmp, "samp1", "storage", "bin_stats_ext.tsv"
        )
        with self.assertRaisesRegex(
            FileNotFoundError, f'file {missing_report_fp} could not be found'
        ):
            _evaluate_bins(
                results_dir=self._tmp, bins=self.bins, db_path=self.db_path,
                common_args=['--reduced_tree', '--threads', '2']
            )


if __name__ == '__main__':
    unittest.main()
