# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import filecmp

import pandas as pd
import pandas.testing as pdt
import qiime2
from qiime2.plugin.testing import TestPluginBase
from qiime2.sdk.parallel_config import (
    NON_QIIMETEST_TEST_CONFIG, ParallelConfig, load_config_from_dict
)

from q2_moshpit.eggnog import (
    _eggnog_annotate
)
from q2_types.genome_data import (OrthologAnnotationDirFmt,
                                  SeedOrthologDirFmt, OrthologFileFmt)
from q2_types.reference_db import (
    EggnogRefDirFmt
)


class TestAnnotate(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'
    parallel_config, mapping = load_config_from_dict(NON_QIIMETEST_TEST_CONFIG)

    def setUp(self):
        super().setUp()
        self.eggnog_db = EggnogRefDirFmt(
            self.get_data_path('eggnog_db/'), mode='r'
        )
        self.eggnog_db_artifact = qiime2.Artifact.import_data(
            'ReferenceDB[Eggnog]',
            self.get_data_path('eggnog_db/')
        )
        self.eggnog_annotate = \
            self.plugin.pipelines["eggnog_annotate"]
        self._eggnog_annotate = \
            self.plugin.methods["_eggnog_annotate"]

    def test_small_good_hits(self):
        seed_orthologs = SeedOrthologDirFmt(
            self.get_data_path('good_hits/'), mode='r'
        )

        obs_obj = _eggnog_annotate(
            eggnog_hits=seed_orthologs, eggnog_db=self.eggnog_db
        )

        exp_fp = self.get_data_path(
            'expected/test_output.emapper.annotations'
        )
        exp = OrthologFileFmt(exp_fp, mode='r').view(pd.DataFrame)

        objs = list(obs_obj.annotations.iter_views(OrthologFileFmt))
        self.assertEqual(len(objs), 1)
        df = objs[0][1].view(pd.DataFrame)
        pdt.assert_frame_equal(df, exp)

    def test_eggnog_annotate_parallel(self):
        orthologs = qiime2.Artifact.import_data(
            'SampleData[Orthologs]',
            self.get_data_path('good_hits/')
        )

        with ParallelConfig(parallel_config=self.parallel_config,
                            action_executor_mapping=self.mapping):
            parallel, = self.eggnog_annotate.parallel(
                    orthologs,
                    self.eggnog_db_artifact
                )._result()

        single, = self._eggnog_annotate(
            eggnog_hits=orthologs,
            eggnog_db=self.eggnog_db_artifact
        )

        parallel = parallel.view(OrthologAnnotationDirFmt)
        single = single.view(OrthologAnnotationDirFmt)

        compare_dir = filecmp.dircmp(parallel.path, single.path)
        self.assertEqual(len(compare_dir.common), 1)

        # TODO: add exact file comparison
