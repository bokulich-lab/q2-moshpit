# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import tempfile
import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.busco.plots_detailed import _draw_detailed_plots
from q2_moshpit.busco.plots_summary import _draw_marker_summary_histograms, \
    _draw_selectable_summary_histograms


class TestBUSCOPlots(TestPluginBase):
    package = "q2_moshpit.busco.tests"

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.df = pd.read_csv(
            self.get_data_path('summaries/all_renamed_with_lengths.csv')
        )

    def test_draw_detailed_plots(self):
        obs = _draw_detailed_plots(
            df=self.df,
            width=100,
            height=250,
            label_font_size=10,
            title_font_size=15,
            spacing=5
        )

        self.assertIsInstance(obs, dict)
        self.assertIn('config', obs)
        self.assertIn('hconcat', obs)
        self.assertIn('data', obs['hconcat'][0])
        for i in (0, 1):
            self.assertEqual(
                obs['hconcat'][i]['spacing'], 5
            )
            self.assertEqual(
                obs['hconcat'][i]['spec']['height']['step'], 250
            )
            self.assertEqual(
                obs['hconcat'][i]['spec']['width'], 100
            )

        config = obs['config']
        for key in ['axis', 'header', 'legend']:
            self.assertEqual(config[key]['labelFontSize'], 10)
            self.assertEqual(config[key]['titleFontSize'], 15)

    def test_draw_marker_summary_histograms(self):
        obs = _draw_marker_summary_histograms(data=self.df)

        self.assertIsInstance(obs, dict)
        self.assertIn('config', obs)
        self.assertIn('vconcat', obs)
        self.assertEqual(len(obs['vconcat'][0]['hconcat']), 4)
        self.assertEqual(len(obs['vconcat'][1]['hconcat']), 4)

        exp_titles = ['Single', 'Duplicated', 'Fragmented', 'Missing']
        obs_titles = [
            x['encoding']['x']['title']
            for x in obs['vconcat'][0]['hconcat']
        ]
        self.assertListEqual(exp_titles, obs_titles)

        exp_titles = ['Scaffolds', 'Contigs n50', 'Scaffold n50', 'Length']
        obs_titles = [
            x['encoding']['x']['title']
            for x in obs['vconcat'][1]['hconcat']
        ]
        self.assertListEqual(exp_titles, obs_titles)

        config = obs['config']
        for key in ['axis', 'header', 'legend']:
            self.assertEqual(config[key]['labelFontSize'], 12)
            self.assertEqual(config[key]['titleFontSize'], 15)

    def test_draw_selectable_summary_histograms(self):
        obs = _draw_selectable_summary_histograms(data=self.df)

        self.assertIsInstance(obs, dict)
        self.assertIn('config', obs)
        self.assertIn('mark', obs)
        self.assertIn('transform', obs)
        self.assertIn('filter', obs['transform'][0])
        self.assertEqual(obs['mark']['type'], 'bar')

        config = obs['config']
        for key in ['axis', 'header', 'legend']:
            self.assertEqual(config[key]['labelFontSize'], 12)
            self.assertEqual(config[key]['titleFontSize'], 15)
