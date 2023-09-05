# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.testing as pdt

import qiime2
from qiime2.plugin.testing import TestPluginBase
from .._method import eggnog_diamond_search, eggnog_annotate
from q2_types_genomics.reference_db import (
    DiamondDatabaseDirFmt, EggnogRefDirFmt)
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from q2_types_genomics.genome_data import SeedOrthologDirFmt, OrthologFileFmt


class TestDiamond(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

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

        pdt.assert_frame_equal(obs, exp)


class TestAnnotate(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    def test_small_good_hits(self):
        so_fp = self.get_data_path('good_hits/')
        seed_orthologs = SeedOrthologDirFmt(so_fp, mode='r')

        egg_db_fp = self.get_data_path('eggnog_db/')
        egg_db = EggnogRefDirFmt(egg_db_fp, mode='r')

        obs_obj = eggnog_annotate(eggnog_hits=seed_orthologs, eggnog_db=egg_db)

        exp_fp = self.get_data_path('expected/test_output.emapper.annotations')
        exp = OrthologFileFmt(exp_fp, mode='r').view(pd.DataFrame)

        objs = list(obs_obj.annotations.iter_views(OrthologFileFmt))
        self.assertEqual(len(objs), 1)
        df = objs[0][1].view(pd.DataFrame)
        pdt.assert_frame_equal(df, exp)
