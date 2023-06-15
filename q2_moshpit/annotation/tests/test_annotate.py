# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from qiime2.plugin.testing import TestPluginBase
from .._method import eggnog_annotate_seed_orthologs
from q2_types_genomics.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from q2_types_genomics.reference_db import EggnogRefDirFmt
import pandas as pd
import pandas.testing as pdt


class TestAnnotate(TestPluginBase):
    package = 'q2_moshpit.annotation.tests'

    def test_small_good_hits(self):
        so_fp = self.get_data_path('good_hits/')
        seed_orthologs = SeedOrthologDirFmt(so_fp, mode='r')

        egg_db_fp = self.get_data_path('eggnog_db/')
        egg_db = EggnogRefDirFmt(egg_db_fp, mode='r')

        obs_obj = eggnog_annotate_seed_orthologs(hits_table=seed_orthologs,
                                                 eggnog_db=egg_db)

        exp_fp = self.get_data_path('expected/test_output.emapper.annotations')
        exp = OrthologFileFmt(exp_fp, mode='r').view(pd.DataFrame)

        objs = list(obs_obj.annotations.iter_views(OrthologFileFmt))
        self.assertEqual(len(objs), 1)
        df = objs[0][1].view(pd.DataFrame)
        pdt.assert_frame_equal(df, exp)
