# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.busco.types import BUSCOResults, BUSCOResultsDirectoryFormat


class TestBUSCOTypes(TestPluginBase):
    package = "q2_moshpit.busco.types.tests"

    def test_feature_data_semantic_type_registration(self):
        self.assertRegisteredSemanticType(BUSCOResults)

    def test_sequence_semantic_type_to_format_registration(self):
        self.assertSemanticTypeRegisteredToFormat(
            BUSCOResults, BUSCOResultsDirectoryFormat
        )
