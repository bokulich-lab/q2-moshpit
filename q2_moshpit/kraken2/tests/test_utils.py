# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest

from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.kraken2.utils import _process_kraken2_arg


class TestKraken2Utils(TestPluginBase):
    package = "q2_moshpit.kraken2.tests"

    def test_process_kraken2_arg_bool(self):
        obs = _process_kraken2_arg("quick", True)
        exp = ["--quick"]
        self.assertListEqual(obs, exp)

    def test_process_kraken2_arg_number(self):
        obs = _process_kraken2_arg("threads", 3)
        exp = ["--threads", "3"]
        self.assertListEqual(obs, exp)

    def test_process_kraken2_arg_string(self):
        obs = _process_kraken2_arg("db", "/some/where/test_db")
        exp = ["--db", "/some/where/test_db"]
        self.assertListEqual(obs, exp)

    def test_process_kraken2_arg_unknown(self):
        with self.assertRaisesRegex(
            NotImplementedError,
            "Parsing arguments of type \"<class 'list'>\" " "is not supported.",
        ):
            _process_kraken2_arg("fake_param", [1, 2])


if __name__ == "__main__":
    unittest.main()
