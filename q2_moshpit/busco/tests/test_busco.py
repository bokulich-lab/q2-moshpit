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
import contextlib
import pandas as pd
from q2_moshpit.busco.busco import evaluate_busco
from q2_moshpit.busco.utils import (
    _parse_busco_params,
    _draw_busco_plots,
    _zip_busco_plots,
    _run_busco,
    _draw_busco_plots_for_render,
    _collect_summaries_and_save,
    _parse_df_columns,
)
from unittest.mock import patch, call, ANY
from qiime2.plugin.testing import TestPluginBase
from q2_types_genomics.per_sample_data._format import MultiMAGSequencesDirFmt


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

    def setUp(self):
        super().setUp()
        with contextlib.ExitStack() as stack:
            self._tmp = stack.enter_context(tempfile.TemporaryDirectory())
            self.addCleanup(stack.pop_all().close)

    # Test `_parse_busco_params`
    def test_parse_busco_params_1(self):
        observed = _parse_busco_params("auto_lineage", True)
        expected = ["--auto-lineage"]
        self.assertSetEqual(set(observed), set(expected))

    def test_parse_busco_params_2(self):
        observed = _parse_busco_params("evalue", 0.66)
        expected = ["--evalue", str(0.66)]
        self.assertSetEqual(set(observed), set(expected))

    def test_parse_busco_params_3(self):
        observed = _parse_busco_params("augustus", True)
        expected = ["--augustus"]
        self.assertSetEqual(set(observed), set(expected))

    def test_parse_busco_params_4(self):
        observed = _parse_busco_params("lineage_dataset", "bacteria-XYZ")
        expected = ["--lineage_dataset", "bacteria-XYZ"]
        self.assertSetEqual(set(observed), set(expected))

    def test_collect_summaries_and_save(self):
        """
        Test for `_collect_summaries_and_save` function.
        Uses data stored in ./data. Checks for data frame equality.
        """
        with tempfile.TemporaryDirectory() as tmp_path:
            path_to_summaries = {}

            for i in range(1, 4):
                path_to_summaries[f"sample{i}"] = self.get_data_path(
                    filename=f"batch_summary_sample{i}.txt"
                )

            observed = _collect_summaries_and_save(
                path_to_run_summaries=path_to_summaries,
                all_summaries_path=os.path.join(tmp_path, "aggregated.csv"),
            )

            expected = pd.read_csv(
                self.get_data_path(filename="all_batch_summaries.csv")
            )
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
        p = self.get_data_path(f"{filename}")
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
        p = self.get_data_path("all_batch_summaries.csv")
        all_summaries_df = pd.read_csv(p)

        # Draw plot
        observed = _draw_busco_plots_for_render(
            all_summaries_df,
            width=600,
            height=18,
            titleFontSize=20,
            labelFontSize=17,
        )

        # Load expected data
        p = self.get_data_path("plot_as_dict.json")
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
        output_dir = self.get_data_path("busco_output")
        sample_ids = os.listdir(output_dir)

        # Initialize assertion objects
        expected = {}
        calls = []

        # Define command arguments
        fake_props = ["--a", "--b", "0.6"]

        # Fabricate list of calls and the expected output
        for sample_id in sample_ids:
            # Make a dictionary to compare output
            p = os.path.join(output_dir, sample_id, "batch_summary.txt")
            expected[sample_id] = p

            # Append call to list of calls to assert the patching
            calls.append(call(
                [
                    "busco",
                    "--a",
                    "--b", "0.6",
                    "--in", self.get_data_path(f"{sample_id}"),
                    "--out_path", output_dir,
                    "-o", sample_id
                ],
                check=True
            ))

        # Run busco and save paths to run summaries
        observed = _run_busco(
            output_dir=output_dir,
            mags=self.mags,
            params=fake_props,
        )

        # Assert output
        self.assertDictEqual(expected, observed)

        # Check for appropiate calls
        subp_run.assert_has_calls(calls, any_order=True)

    @patch("subprocess.run")
    def test_run_busco_exception(self, subp_run):
        """
        Test function `_run_busco`. Checks for a raised exception.
        """
        with tempfile.TemporaryDirectory() as tmp_path:
            # Define command arguments
            fake_props = ["--a", "--b", "0.6"]
            output_dir = os.path.join(tmp_path, "busco_output")

            with self.assertRaises(FileNotFoundError):
                # Run busco and save paths to run summaries
                _ = _run_busco(
                    output_dir=output_dir,
                    mags=self.mags,
                    params=fake_props,
                )

        # Assert that the patch was called once.
        cmd = [
            "busco",
            "--a",
            "--b", "0.6",
            "--in", self.get_data_path("sample1"),
            "--out_path", output_dir,
            "-o", "sample1"
        ]
        subp_run.assert_called_once_with(cmd, check=True)

    # Integration test busco.
    @patch('q2_moshpit.busco.utils._run_busco')
    @patch('q2_moshpit.busco.utils._zip_busco_plots')
    @patch('q2_moshpit.busco.utils._draw_busco_plots')
    @patch('q2_moshpit.busco.utils._collect_summaries_and_save')
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
                `_collect_summaries_and_save`
            zip_busco_plots (unittest.mock): mock object for function
                `_draw_busco_plots`.
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

    def test_parse_df_columns(self):
        # This side effect will return the all_summaries_dfs
        p1 = self.get_data_path("all_batch_summaries.csv")
        observed = pd.read_csv(p1)
        observed = _parse_df_columns(observed)

        p2 = self.get_data_path("all_batch_summaries_formatted.csv")
        expected = pd.read_csv(p2)

        try:
            pd.testing.assert_frame_equal(observed, expected)
        except AssertionError as e:
            print(e)
            self.assertTrue(False)
        else:
            self.assertTrue(True)
