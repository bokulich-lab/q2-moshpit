"""
Unit test for BUSCO wrapper.
"""
import os
from ..busco import _parse_busco_params


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


# TODO: Test `_zip_busco_plots`
# TODO: Test `_draw_busco_plots`
# TODO: Test `_run_busco`
# TODO: Test `_draw_busco_plots_for_render`
# TODO: Integration test busco.
