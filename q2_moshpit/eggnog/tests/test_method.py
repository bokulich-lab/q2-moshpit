# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import filecmp
import qiime2
import pandas as pd
import pandas.testing as pdt
from qiime2.plugin.testing import TestPluginBase
from qiime2.sdk.parallel_config import ParallelConfig
from q2_types.feature_data_mag import MAGSequencesDirFmt
from .._method import _eggnog_diamond_search, _eggnog_annotate
from q2_types.reference_db import (
    DiamondDatabaseDirFmt, EggnogRefDirFmt
)
from q2_types.per_sample_sequences import ContigSequencesDirFmt
from q2_types.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from q2_types.feature_data_mag import OrthologAnnotationDirFmt


class TestDiamond(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    def setUp(self):
        super().setUp()
        self.diamond_db_artifact = qiime2.Artifact.import_data(
            'ReferenceDB[Diamond]',
            self.get_data_path('random-db-1')
        )
        self.diamond_db = self.diamond_db_artifact.view(DiamondDatabaseDirFmt)

        self.eggnog_diamond_search = \
            self.plugin.pipelines["eggnog_diamond_search"]
        self._eggnog_diamond_search = \
            self.plugin.methods["_eggnog_diamond_search"]

    def test_good_small_search_contigs(self):
        contigs = qiime2.Artifact.import_data(
            'SampleData[Contigs]',
            self.get_data_path('contig-sequences-1')
        ).view(ContigSequencesDirFmt)

        _, obs = _eggnog_diamond_search(
            sequences=contigs,
            diamond_db=self.diamond_db
        )
        exp = pd.DataFrame({'0': [1.0, 0.0], '2': [1.0, 0.0], '8': [0.0, 3.0]},
                           index=['s1', 's2'])
        exp.columns.name = 'sseqid'

        pdt.assert_frame_equal(obs, exp)

    def test_good_small_search_mags(self):
        mags = qiime2.Artifact.import_data(
            'FeatureData[MAG]',
            self.get_data_path('mag-sequences')
        ).view(MAGSequencesDirFmt)

        _, obs = _eggnog_diamond_search(
            sequences=mags,
            diamond_db=self.diamond_db
        )
        exp = pd.DataFrame(
            {'8': [3.0, 0.0], '0': [0.0, 1.0], '2': [0.0, 1.0]},
            index=['194b1aca-9373-4298-ba5c-05b382d1f553',
                   'e1ba1d20-b466-4fef-ae19-6bf9c5c63d6f']
        )
        exp.columns.name = 'sseqid'

        pdt.assert_frame_equal(obs, exp)

    def test_eggnog_search_parallel_contigs(self):
        contigs = qiime2.Artifact.import_data(
            'SampleData[Contigs]',
            self.get_data_path('contig-sequences-1')
        )

        with ParallelConfig():
            _, parallel = self.eggnog_diamond_search.parallel(
                    contigs,
                    self.diamond_db_artifact
                )._result()

        _, single = self._eggnog_diamond_search(
            sequences=contigs,
            diamond_db=self.diamond_db_artifact
        )

        parallel = parallel.view(pd.DataFrame)
        single = single.view(pd.DataFrame)

        pdt.assert_frame_equal(parallel, single)

    def test_eggnog_search_parallel_mags(self):
        mags = qiime2.Artifact.import_data(
            'FeatureData[MAG]',
            self.get_data_path('mag-sequences')
        )

        with ParallelConfig():
            _, parallel = self.eggnog_diamond_search.parallel(
                    mags,
                    self.diamond_db_artifact
                )._result()

        _, single = self._eggnog_diamond_search(
            sequences=mags,
            diamond_db=self.diamond_db_artifact
        )

        parallel = parallel.view(pd.DataFrame)
        single = single.view(pd.DataFrame)

        pdt.assert_frame_equal(parallel, single)


class TestAnnotate(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

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
            'SampleData[BLAST6]',
            self.get_data_path('good_hits/')
        )

        with ParallelConfig():
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
