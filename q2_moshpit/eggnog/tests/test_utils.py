# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_moshpit.eggnog._utils import _check_taxa
from qiime2.plugin.testing import TestPluginBase


class TestDButils(TestPluginBase):
    package = 'q2_moshpit.eggnog'

    def test_raise_on_mixed_numeric(self):
        mixed_taxa = ["1184", "38", "Streptophyta"]
        with self.assertRaisesRegex(ValueError,
                                    "All taxa inputs must be the same type,"
                                    " either all taxids as integers or all"
                                    " string labels."):
            _check_taxa(mixed_taxa)

    def test_raise_on_mixed_string(self):
        mixed_taxa = ["Streptophyta", "1184", "38"]
        with self.assertRaisesRegex(ValueError,
                                    "All taxa inputs must be the same type,"
                                    " either all taxids as integers or all"
                                    " string labels."):
            _check_taxa(mixed_taxa)

    def test_successful_numeric(self):
        # smallest taxa download available
        good_numeric = ["85004"]
        _check_taxa(good_numeric)

    def test_successful_string(self):
        good_str = ["Bifidobacteriales"]
        _check_taxa(good_str)
