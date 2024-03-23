# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import pandas as pd
from q2_moshpit.busco.busco import evaluate_busco
from unittest.mock import patch, ANY
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences._format import MultiMAGSequencesDirFmt


class TestBUSCO(TestPluginBase):
    package = "q2_moshpit.busco.tests"

    @classmethod
    def setUpClass(self):

        # Set base path
        p = os.path.join(os.path.dirname(__file__), "data")

        # Get MAGs fixture
        self.mags = MultiMAGSequencesDirFmt(
            path=p,
            mode="r",
        )

    # Integration test busco.
    @patch('q2_moshpit.busco.utils._run_busco')
    @patch('q2_moshpit.busco.utils._zip_busco_plots')
    @patch('q2_moshpit.busco.utils._draw_detailed_plots')
    @patch('q2_moshpit.busco.utils._collect_summaries')
    def test_integration_busco(
        self,
        collect_summaries,
        draw_busco_plots,
        zip_busco_plots,
        run_busco
    ):
        """
        Tests entire busco run and patches the previously tested functions.
        Checks for existence of HTML output.

        Args:
            collect_summaries (unittest.mock): mock object for function
                `_collect_summaries`
            zip_busco_plots (unittest.mock): mock object for function
                `_draw_detailed_plots`.
            zip_busco_plots (unittest.mock): mock object for function
                `_zip_busco_plots`.
            run_busco (unittest.mock): mock object for function
                `_run_busco`.
        """
        # import shutil
        # path_to_look_at_html = "/Users/santiago/Downloads/busco_debug_bench"

        with tempfile.TemporaryDirectory() as tmp_path:
            # This side effect will return the all_summaries_dfs
            p = self.get_data_path("all_batch_summaries.csv")
            collect_summaries.return_value = pd.read_csv(p)

            # Run busco
            evaluate_busco(output_dir=str(tmp_path), bins=self.mags)

            # For render debugging
            # shutil.copytree(str(tmp_path), path_to_look_at_html)

            # Check for the existence of the html file
            self.assertTrue(os.path.exists(f"{tmp_path}/index.html"))

            # Assert that the calls where done properly
            run_busco.assert_called_once_with(
                output_dir=run_busco.call_args.kwargs['output_dir'],
                mags=self.mags,
                params=[
                    '--mode', 'genome',
                    '--cpu', '1',
                    '--contig_break', '10',
                    '--evalue', '0.001',
                    '--limit', '3'
                ]
            )

            collect_summaries.assert_called_once_with(
                all_summaries_path=os.path.join(
                    tmp_path, "all_batch_summaries.csv"
                ),
                path_to_run_summaries=ANY,
            )

            draw_busco_plots.assert_called_once_with(
                path_to_run_summaries=ANY,
                plots_dir=draw_busco_plots.call_args.kwargs["plots_dir"]
            )

            zip_busco_plots.assert_called_once_with(
                paths_to_plots=ANY,
                zip_path=os.path.join(tmp_path, "busco_plots.zip")
            )
