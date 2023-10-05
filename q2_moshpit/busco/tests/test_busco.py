# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import zipfile
import pandas as pd
from q2_moshpit.busco.busco import busco
from q2_moshpit.busco.utils import (
    _parse_busco_params,
    _draw_busco_plots,
    _zip_busco_plots,
    _run_busco,
    _draw_busco_plots_for_render,
    _collect_summaries_and_save,
    _parse_df_columns,
)
from unittest.mock import patch
from qiime2.plugin.testing import TestPluginBase
from q2_types_genomics.per_sample_data._format import MultiMAGSequencesDirFmt


class TestBUSCO(TestPluginBase):
    package = "q2_moshpit.busco.tests"

    @classmethod
    def setUpClass(cls):

        # Set base path
        p = os.path.join(os.path.dirname(__file__), "data")

        # Get MAGs fixture
        cls.mags = MultiMAGSequencesDirFmt(
            path=p,
            mode="r",
        )

    # Test `_parse_busco_params`
    def test_parse_busco_params_1(self):
        observed = _parse_busco_params("auto_lineage", True)
        expected = ["--auto-lineage"]
        self.assertSetEqual(set(observed), set(expected))

    def test_parse_busco_params_2(self):
        observed = _parse_busco_params("evalue", 0.66)
        expected = ["--evalue", str(0.66)]
        self.assertSetEqual(set(observed), set(expected))

    def test_collect_summaries_and_save(self):
        """
        Test for `_collect_summaries_and_save` function.
        Uses data stored in ./data. Checks for data frame equality.
        """
        with tempfile.TemporaryDirectory() as tmp_path:
            path_to_summaries = {}
            common_path = os.path.join(os.path.dirname(__file__), "data")

            for i in range(1, 4):
                path_to_summaries[f"sample{i}"] = os.path.join(
                    common_path,
                    f"batch_summary_sample{i}.txt"
                )

            observed = _collect_summaries_and_save(
                path_to_run_summaries=path_to_summaries,
                all_summaries_path=os.path.join(tmp_path, "aggregated.csv"),
            )

            p = os.path.join(common_path, "all_batch_summaries.csv")
            expected = pd.read_csv(p)
            pd.set_option('display.max_columns', None)

            try:
                pd.testing.assert_frame_equal(observed, expected)
            except AssertionError as e:
                print(e)
                self.assertTrue(False)
            else:
                self.assertTrue(True)

    # Test `_draw_busco_plots`
    def draw_n_busco_plots(self, filename, delim):
        """
        Creates plot from a table containing information about
        one or more samples. Checks for the existence of the output
        plots.

        Args:
            filename (str): name of file in ./data to construct the images
            delim (str): delimiter of `filename`
        """
        # Create an empty dictionary to store the DataFrames
        path_to_run_summaries = {}

        # Group the DataFrame by the 'sample_id' column
        p = os.path.join(
            os.path.dirname(__file__),
            f"data/{filename}"
        )
        df = pd.read_csv(p, delimiter=delim)
        grouped = df.groupby("sample_id")

        # Iterate through the groups and store each group as a
        # DataFrame in the dictionary
        # Creates output directory with path 'tmp'
        with tempfile.TemporaryDirectory() as tmp_path:
            for sample_id, group_df in grouped:
                path_to_df = f"{tmp_path}/{sample_id}.csv"
                group_df.to_csv(path_to_df, sep="\t", index=False)
                path_to_run_summaries[sample_id] = path_to_df

            # Draw plots
            paths_to_plots = _draw_busco_plots(
                path_to_run_summaries=path_to_run_summaries,
                plots_dir=os.path.join(tmp_path, "plots"),
            )

            # Check if busco plots are in fact generated
            for _, value in paths_to_plots.items():
                self.assertTrue(os.path.exists(value))

    def test_draw_busco_plots_multiple(self):
        self.draw_n_busco_plots(
            filename="all_batch_summaries.csv", delim=","
        )

    def test_draw_busco_plots_one(self):
        self.draw_n_busco_plots(
            filename="batch_summary_sample1.txt", delim="\t"
        )

    # Test `_draw_busco_plots_for_render`
    def test_draw_busco_plots_for_render(self):
        """
        Tests function `_draw_busco_plots_for_render`.
        Checks for dictionary equality.
        """
        # Load data
        p = os.path.dirname(os.path.abspath(__file__))
        p2 = os.path.join(p, "data/all_batch_summaries.csv")
        all_summaries_df = pd.read_csv(p2)

        # Draw plot
        observed = _draw_busco_plots_for_render(
            all_summaries_df,
            width=600,
            height=18,
            titleFontSize=20,
            labelFontSize=17,
        )

        # Load expected data
        p = os.path.join(p, "data/plot_as_dict.json")
        with open(p, "r") as json_file:
            expected = json_file.read()

        # self.maxDiff = None
        self.assertEqual(expected, observed)

    # Test `_draw_busco_plots`
    def mock_draw_busco_plots(self, tmp_path: str, num_files: int) -> dict:
        """
        Mocks the generation of sample wise plots by generating
        empty files.

        Args:
            tmp_path (str): Path where to write the empty files.
            num_files (int): number of empty files to create, one per sample.

        Returns:
            paths_to_plots (dict):  dictionary with keys sample_id and value
                path to empty file.
        """
        # Generate random images
        paths_to_plots = {}

        # Path to output
        out_dir = os.path.join(tmp_path, "zip_this_dir")
        os.makedirs(out_dir)

        # Loop to create the empty files
        for i in range(num_files):
            # Specify the name of each empty file
            file_name = f"empty_file_{i}.svg"

            # Combine the directory path and file name to create the full
            # file path
            file_path = os.path.join(out_dir, file_name)
            paths_to_plots[f"empty_file_{i}"] = file_path

            # Create an empty file
            with open(file_path, 'w'):
                pass

        return paths_to_plots

    def test_zip_busco_plots_multiple(self):
        """
        Checks for existence of zip file.
        """
        with tempfile.TemporaryDirectory() as tmp_path:
            paths_to_plots = self.mock_draw_busco_plots(
                num_files=6, tmp_path=tmp_path
            )

            # Zip graphs for user download
            zip_path = os.path.join(tmp_path, "busco_plots.zip")
            _zip_busco_plots(paths_to_plots=paths_to_plots, zip_path=zip_path)

            # Check for existence of file
            self.assertTrue(zipfile.is_zipfile(zip_path))

    def test_zip_busco_plots_one(self):
        """
        Checks for existence of zip file.
        """
        with tempfile.TemporaryDirectory() as tmp_path:
            paths_to_plots = self.mock_draw_busco_plots(
                num_files=1, tmp_path=tmp_path
            )

            # Zip graphs for user download
            zip_path = os.path.join(tmp_path, "busco_plots.zip")
            _zip_busco_plots(paths_to_plots=paths_to_plots, zip_path=zip_path)

            # Check for existence of file
            self.assertTrue(zipfile.is_zipfile(zip_path))

    @patch('subprocess.run')
    def test_run_busco(self, subp_run):
        """
        Test function `_run_busco`. Checks for dictionary equality.
        """
        p = os.path.dirname(os.path.abspath(__file__))
        p2 = os.path.join(p, "data", "busco_output")

        sample_ids = os.listdir(p2)

        expected = {}
        for sample_id in sample_ids:
            expected[sample_id] = os.path.join(
                p2, sample_id, "batch_summary.txt"
            )

        # Define command arguments
        fake_props = {"a": "b", "c": "d"}

        # Run busco and save paths to run summaries
        observed = _run_busco(
            output_dir=p2,
            mags=self.mags,
            params=fake_props,
        )

        self.assertDictEqual(expected, observed)

    @patch("subprocess.run")
    def test_run_busco_exception(self, subp_run):
        """
        Test function `_run_busco`. Checks for a raised exception.
        """
        with tempfile.TemporaryDirectory() as tmp_path:
            # Define command arguments
            fake_props = {"a": "b", "c": "d"}

            with self.assertRaises(FileNotFoundError):
                # Run busco and save paths to run summaries
                _ = _run_busco(
                    output_dir=os.path.join(tmp_path, "busco_output"),
                    mags=self.mags,
                    params=fake_props,
                )

    # Integration test busco.
    @patch('q2_moshpit.busco.utils._run_busco')
    @patch('q2_moshpit.busco.utils._zip_busco_plots')
    @patch('q2_moshpit.busco.utils._draw_busco_plots')
    @patch('q2_moshpit.busco.utils._collect_summaries_and_save')
    def test_integration_busco(
        self,
        collect_summaries,
        not_used_1,
        not_used_2,
        not_used_3
    ):
        """
        Tests entire busco run and patches the previously tested functions.
        Checks for existence of HTML output.

        Args:
            collect_summaries (unittest.mock): mock object for function
                `_collect_summaries_and_save`
            not_used_1 (unittest.mock): mock object for function
                `_draw_busco_plots`. Not used.
            not_used_2 (unittest.mock): mock object for function
                `_zip_busco_plots`. Not used.
            not_used_3 (unittest.mock): mock object for function
                `_run_busco`.
        """
        # import shutil
        # path_to_look_at_html = "/Users/santiago/Downloads/busco_debug_bench"

        with tempfile.TemporaryDirectory() as tmp_path:
            # Define side effects and return arguments for patches
            # This side effect will return the all_summaries_dfs
            p = os.path.join(
                os.path.dirname(__file__),
                "data/all_batch_summaries.csv"
            )
            collect_summaries.return_value = pd.read_csv(p)

            # Run busco
            busco(output_dir=str(tmp_path), bins=self.mags)

            # For render debugging
            # shutil.copytree(str(tmp_path), path_to_look_at_html)

            # Check for the existence of the html file
            self.assertTrue(os.path.exists(f"{tmp_path}/index.html"))

    def test_parse_df_columns(self):
        # This side effect will return the all_summaries_dfs
        p1 = os.path.join(
            os.path.dirname(__file__),
            "data/all_batch_summaries.csv"
        )
        p2 = os.path.join(
            os.path.dirname(__file__),
            "data/all_batch_summaries_formatted.csv"
        )
        observed = pd.read_csv(p1)
        observed = _parse_df_columns(observed)
        expected = pd.read_csv(p2)

        try:
            pd.testing.assert_frame_equal(observed, expected)
        except AssertionError as e:
            print(e)
            self.assertTrue(False)
        else:
            self.assertTrue(True)
