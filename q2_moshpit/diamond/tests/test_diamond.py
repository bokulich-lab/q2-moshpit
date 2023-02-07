# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from qiime2.plugin.testing import TestPluginBase
from .._method import eggnog_diamond_search
from q2_types_genomics.reference_db import DiamondDatabaseDirFmt
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt


class TestDiamond(TestPluginBase):
    package = 'q2_moshpit.diamond.tests'

    def test_good_small_search(self):
        input_sequences = ContigSequencesDirFmt(
                self.get_data_path('tiny_test_data.qza'), mode='r')

        diamond_db = DiamondDatabaseDirFmt(
                self.get_data_path('tiny_diamond_db.qza'), mode='r')

        eggnog_diamond_search(
                input_sequences=input_sequences,
                diamond_db=diamond_db)
