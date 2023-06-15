# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd

import qiime2
from qiime2.plugin.testing import TestPluginBase
from .._method import eggnog_diamond_search
from q2_types_genomics.reference_db import DiamondDatabaseDirFmt
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt


class TestDiamond(TestPluginBase):
    package = 'q2_moshpit.diamond.tests'

    def test_good_small_search(self):
        input_sequences = qiime2.Artifact.import_data(
            'SampleData[Contigs]',
            self.get_data_path('contig-sequences-1')
        ).view(ContigSequencesDirFmt)

        diamond_db = qiime2.Artifact.import_data(
            'ReferenceDB[Diamond]',
            self.get_data_path('random-db-1')
        ).view(DiamondDatabaseDirFmt)

        _, obs = eggnog_diamond_search(
                input_sequences=input_sequences,
                diamond_db=diamond_db)
        exp = pd.DataFrame({'0': [1.0, 0.0], '2': [1.0, 0.0], '8': [0.0, 3.0]},
                           index=['s1', 's2'])
        exp.columns.name = 'sseqid'

        pd.testing.assert_frame_equal(obs, exp)
