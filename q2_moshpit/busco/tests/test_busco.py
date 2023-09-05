"""
Unit test for BUSCO wrapper.
"""
import os
import json
import pandas as pd
from ..busco import (
    _parse_busco_params,
    _draw_busco_plots,
    _zip_busco_plots,
    _run_busco,
    _draw_busco_plots_for_render,
    busco,
)
from ..._utils import _process_common_input_params


# Test that fixtrue is working.
def test_tmp_path_factory_fixtrue(MAGS_test_data):
    assert os.path.exists(MAGS_test_data.path)


# Test `_parse_busco_params`
def test_parse_busco_params_1():
    observed = _parse_busco_params("auto_lineage", True)
    expected = ["--auto-lineage"]
    assert set(observed) == set(expected)


def test_parse_busco_params_2():
    observed = _parse_busco_params("evalue", 0.66)
    expected = ["--evalue", str(0.66)]
    assert set(observed) == set(expected)


def test_process_common_input_params():
    data = {
        "a": 0,
        "b": 1,
        "c": 0.0,
        "d": 3.14,
        "e": "",
        "f": "Hello",
        "g": None,
        "h": 42,
        "i": 0.0,
        "j": [],
        "k": "World",
        "l": False,
        "m": True,
    }
    observed = _process_common_input_params(_parse_busco_params, data)
    expected = [
        "--a",
        "0",
        "--b",
        "1",
        "--c",
        "0.0",
        "--d",
        "3.14",
        "--f",
        "Hello",
        "--h",
        "42",
        "--i",
        "0.0",
        "--k",
        "World",
        "--m",
    ]
    assert set(observed) == set(expected)


# Test `_zip_busco_plots` and `_draw_busco_plots`
def test_zip_busco_plots(tmp_path):
    # Create an empty dictionary to store the DataFrames
    path_to_run_summeries = {}

    # Group the DataFrame by the 'sample_id' column
    p = os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(p, "./data/all_batch_summeries.csv")
    df = pd.read_csv(p)
    grouped = df.groupby("sample_id")

    # Iterate through the groups and store each group as a
    # DataFrame in the dictionary
    for sample_id, group_df in grouped:
        path_to_df = f"{tmp_path}/{sample_id}.csv"
        group_df.to_csv(path_to_df, sep="\t", index=False)
        path_to_run_summeries[sample_id] = path_to_df

    # Draw plots
    paths_to_plots = _draw_busco_plots(
        path_to_run_summeries=path_to_run_summeries,
        output_dir=os.path.join(tmp_path, "plots"),
    )

    # Zip graphs for user download
    zip_name = os.path.join(tmp_path, "busco_plots.zip")
    _zip_busco_plots(paths_to_plots=paths_to_plots, zip_path=zip_name)

    # Check for existence of file
    assert os.path.exists(zip_name)


# Test `_run_busco`
def test_run_busco(tmp_path, MAGS_test_data):
    # Define comand arguments and parse them
    kwargs = {
        "mode": "genome",
        "lineage_dataset": "bacteria_odb10",
        "download_path": tmp_path,
    }
    common_args = _process_common_input_params(
        processing_func=_parse_busco_params, params=kwargs
    )

    # Run busco and save paths to run summeries
    path_to_run_summeries = _run_busco(
        output_dir=os.path.join(tmp_path, "busco_output"),
        mags=MAGS_test_data,
        params=common_args,
    )

    # Test existence of output files
    for _, path_to_summary in path_to_run_summeries.items():
        assert os.path.exists(path_to_summary)


# Test `_draw_busco_plots_for_render`
def test_draw_busco_plots_for_render():
    p = os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(p, "./data/all_batch_summeries.csv")
    all_summeries_df = pd.read_csv(p)
    plot_as_dict = _draw_busco_plots_for_render(
        all_summeries_df,
        width=600,
        height=18,
        titleFontSize=20,
        labelFontSize=17,
    )
    assert isinstance(json.loads(plot_as_dict), dict)


# Integration test busco.
def test_integration_busco(tmp_path, MAGS_test_data):
    # Defin kwargs
    kwargs = {
        "mode": "genome",
        "lineage_dataset": "bacteria_odb10",
        "download_path": tmp_path,
        "output_dir": str(tmp_path),  # Needs to be string
    }

    # Run busco
    busco(**kwargs, bins=MAGS_test_data)

    # Check for the existence of the html file
    assert os.path.exists(f"{tmp_path}/index.html")
