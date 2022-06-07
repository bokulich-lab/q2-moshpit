# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest

import altair as alt
import pandas as pd
from altair import (
    Color, X, Y, FilterTransform, Undefined, Scale, Tooltip, MarkDef
)
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.checkm.plots import (_prep_bar_plot, _prep_scatter_plot,
                                     _prep_contig_plots)


class TestCheckMPlots(TestPluginBase):
    package = 'q2_moshpit.checkm.tests'

    def setUp(self):
        super().setUp()
        self.data = pd.DataFrame({
            'samples': ['s1', 's2'], 'bins': ['b1', 'b1'],
            'metric1': [0.1, 0.5], 'metric2': [58.0, 99.5]
        })
        self.base_plot = alt.Chart(self.data)
        self.filter = alt.selection_single(
            fields=['samples'],
            bind=alt.binding_select(options=['s1', 's2'], name='Sample')
        )
        self.selection = alt.selection_interval()

    def assertPlot(self, obs_spec, exp_spec):
        color = exp_spec.get('color')
        if color:
            self.assertEqual(exp_spec['color'], obs_spec.encoding.color)
        self.assertEqual(exp_spec['x'], obs_spec.encoding.x)
        self.assertEqual(exp_spec['y'], obs_spec.encoding.y)
        self.assertEqual(exp_spec['mark'], obs_spec.mark)
        if exp_spec.get('selection_count'):
            self.assertEqual(
                len(obs_spec.selection), exp_spec['selection_count']
            )
        else:
            self.assertEqual(obs_spec.selection, Undefined)
        if exp_spec.get('transform'):
            self.assertEqual(len(obs_spec.transform), exp_spec['filter_count'])
            for trans in obs_spec.transform:
                self.assertIsInstance(trans, exp_spec['transform'])
        else:
            self.assertEqual(obs_spec.transform, Undefined)
        tooltips = exp_spec.get('tooltip')
        if tooltips:
            [self.assertIn(x, obs_spec.encoding.tooltips) for x in tooltips]

    def test_prep_bar_plot(self):
        obs_spec = _prep_bar_plot(
            self.base_plot, self.filter, x_col='bins', y_col='metric1',
            x_title='Bin', y_title='Metric1', color_shorthand='metric1:Q',
            color_title='Metric1', color_map='blues'
        )
        exp_spec = {
            'color': Color(
                'metric1:Q', title='Metric1',
                scale=alt.Scale(scheme='blues'),
                sort='ascending'
            ),
            'x': X('bins', title='Bin'),
            'y': Y('metric1', title='Metric1'),
            'mark': 'bar',
            'transform': FilterTransform,
            'filter_count': 1,
            'selection': None,
            'selection_count': None
        }
        self.assertPlot(obs_spec, exp_spec)

    def test_prep_bar_plot_with_selection(self):
        obs_spec = _prep_bar_plot(
            self.base_plot, self.filter, x_col='bins', y_col='metric1',
            x_title='Bin', y_title='Metric1', color_shorthand='metric1:Q',
            color_title='Metric1', color_map='blues',
            bin_selection=self.selection
        )
        exp_spec = {
            'color': Color(
                'metric1:Q', title='Metric1',
                scale=alt.Scale(scheme='blues'),
                sort='ascending'
            ),
            'x': X('bins', title='Bin'),
            'y': Y('metric1', title='Metric1'),
            'mark': 'bar',
            'transform': FilterTransform,
            'filter_count': 2,
            'selection_count': None,
        }
        self.assertPlot(obs_spec, exp_spec)

    def test_prep_scatter_plot(self):
        obs_spec = _prep_scatter_plot(
            self.base_plot, x_col='points', y_col='metric1', x_title='Points',
            y_title='Metric1', primary_selection=self.selection,
            selection_col='bin:N', selection_title='Bin', color_map='blues',
            reverse_x=True
        )
        exp_spec = {
            'color': None,
            'x': X('points:Q', title='Points', scale=Scale(reverse=True)),
            'y': Y('metric1:Q', title='Metric1', scale=Scale(reverse=False)),
            'mark': MarkDef(
                filled=True, opacity=0.7, size=80, stroke='lightgrey',
                strokeOpacity=0.8, strokeWidth=0.2, type='point'
            ),
            'transform': None,
            'filter_count': None,
            'selection_count': 2,
            'tooltips': [
                Tooltip('points:Q', title='Points', format='.2'),
                Tooltip('metric1:Q', title='Metric1', format='.2')
            ]
        }
        self.assertPlot(obs_spec, exp_spec)

    def test_prep_scatter_plot_with_filter(self):
        obs_spec = _prep_scatter_plot(
            self.base_plot, x_col='points', y_col='metric1', x_title='Points',
            y_title='Metric1', primary_selection=self.selection,
            selection_col='bin:N', selection_title='Bin', color_map='blues',
            primary_filter=self.filter
        )
        exp_spec = {
            'color': None,
            'x': X('points:Q', title='Points', scale=Scale(reverse=False)),
            'y': Y('metric1:Q', title='Metric1', scale=Scale(reverse=False)),
            'mark': MarkDef(
                filled=True, opacity=0.7, size=80, stroke='lightgrey',
                strokeOpacity=0.8, strokeWidth=0.2, type='point'
            ),
            'transform': FilterTransform,
            'filter_count': 1,
            'selection_count': 3,
            'tooltips': [
                Tooltip('points:Q', title='Points', format='.2'),
                Tooltip('metric1:Q', title='Metric1', format='.2')
            ]
        }
        self.assertPlot(obs_spec, exp_spec)

    def test_prep_contig_plots(self):
        obs_specs = _prep_contig_plots(
            self.base_plot, self.selection, self.filter
        )
        exp_cols = {
            'longest_contig': 'Longest contig length [bp]',
            'n50_contigs': 'N50 contigs [bp]',
            'mean_contig_length': 'Mean contig length [bp]',
            'ambiguous_bases': 'Count of ambiguous bases'
        }
        exp_specs = [{
            'color': None,
            'mark': 'bar',
            'x': X('bin_id:N', title='Bin'),
            'y': Y(f'{k}:Q', title=v),
            'transform': FilterTransform,
            'filter_count': 2,
            'selection_count': None,
        } for k, v in exp_cols.items()]
        for obs_spec, exp_spec in zip(obs_specs.values(), exp_specs):
            self.assertPlot(obs_spec, exp_spec)


if __name__ == '__main__':
    unittest.main()
