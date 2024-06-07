# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data_mag import MAGSequencesDirFmt
from .._utils import (
    _construct_param, _process_common_input_params,
    _calculate_md5_from_file, get_feature_lengths
)


def fake_processing_func(key, val):
    if isinstance(val, bool):
        return [_construct_param(key)]
    else:
        return [_construct_param(key), str(val)]


def fake_processing_func_no_falsy_filtering(key, val):
    """
    NOTE: There is a need for a function that does this since
    `_process_common_input_params` already filter falsy values.
    If a second filter is applied then some parameters are omitted.
    """
    if isinstance(val, bool):
        return [_construct_param(key)]
    else:
        return [_construct_param(key), str(val)]


class TestUtils(TestPluginBase):
    package = 'q2_moshpit.tests'

    def test_construct_param_simple(self):
        obs = _construct_param('test')
        exp = '--test'
        self.assertEqual(obs, exp)

    def test_construct_param_complex(self):
        obs = _construct_param('test_param')
        exp = '--test-param'
        self.assertEqual(obs, exp)

    def test_process_common_inputs_bools(self):
        kwargs = {'arg1': False, 'arg2': True}
        obs = _process_common_input_params(fake_processing_func, kwargs)
        exp = ['--arg2']
        self.assertListEqual(obs, exp)

    def test_process_common_inputs_nones(self):
        kwargs = {'arg1': 'some-value', 'arg2': None}
        obs = _process_common_input_params(fake_processing_func, kwargs)
        exp = ['--arg1', 'some-value']
        self.assertListEqual(obs, exp)

    def test_process_common_inputs_with_values(self):
        kwargs = {'arg1': 'value1', 'arg2': 'value2'}
        obs = _process_common_input_params(fake_processing_func, kwargs)
        exp = ['--arg1', 'value1', '--arg2', 'value2']
        self.assertListEqual(obs, exp)

    def test_process_common_inputs_mix(self):
        kwargs = {
            'arg1': None, 'arg2': 'some-value', 'arg3': False, 'arg4': True
        }
        obs = _process_common_input_params(fake_processing_func, kwargs)
        exp = ['--arg2', 'some-value', '--arg4']
        self.assertListEqual(obs, exp)

    def test_process_common_inputs_mix_with_falsy_values(self):
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
        observed = _process_common_input_params(
            fake_processing_func_no_falsy_filtering, data
        )
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
        self.assertSetEqual(set(observed), set(expected))

    def test_calculate_md5_from_pass(self):
        path_to_file = self.get_data_path("md5/a.txt")
        observed_hash = _calculate_md5_from_file(path_to_file)
        self.assertEqual(observed_hash, "a583054a9831a6e7cc56ea5cd9cac40a")

    def test_calculate_md5_from_fail(self):
        path_to_file = self.get_data_path("md5/b.txt")
        observed_hash = _calculate_md5_from_file(path_to_file)
        self.assertNotEqual(observed_hash, "a583054a9831a6e7cc56ea5cd9cac40a")

    def test_get_feature_lengths(self):
        mags = MAGSequencesDirFmt(self.get_data_path('mags-derep'), mode='r')
        obs = get_feature_lengths(mags)
        exp = pd.DataFrame({
            'id': [
                '24dee6fe-9b84-45bb-8145-de7b092533a1',
                'ca7012fc-ba65-40c3-84f5-05aa478a7585',
                'd65a71fa-4279-4588-b937-0747ed5d604d'
            ],
            'length': [66, 70, 363]
        })
        exp.set_index('id', inplace=True)
        pd.testing.assert_frame_equal(obs, exp)
