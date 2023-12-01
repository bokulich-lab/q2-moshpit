# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import pandas.testing as pdt
import os
from unittest.mock import patch
import qiime2
from qiime2.plugin.testing import TestPluginBase
from q2_types_genomics.feature_data import MAGSequencesDirFmt
from .._method import (
  eggnog_diamond_search, eggnog_annotate, fetch_eggnog_db, build_diamond_db
)
from q2_types_genomics.reference_db import (
    DiamondDatabaseDirFmt, EggnogRefDirFmt)
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from q2_types_genomics.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from q2_types.feature_data import ProteinSequencesDirectoryFormat


class TestDiamond(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    def setUp(self):
        super().setUp()
        self.diamond_db = qiime2.Artifact.import_data(
            'ReferenceDB[Diamond]',
            self.get_data_path('random-db-1')
        ).view(DiamondDatabaseDirFmt)

    def test_good_small_search_contigs(self):
        contigs = qiime2.Artifact.import_data(
            'SampleData[Contigs]',
            self.get_data_path('contig-sequences-1')
        ).view(ContigSequencesDirFmt)

        _, obs = eggnog_diamond_search(
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

        _, obs = eggnog_diamond_search(
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


class TestBuildDiamondDB(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    @patch("subprocess.run")
    def test_build_diamond_db(self, subp_run):
        # Instantiate input
        sequences = ProteinSequencesDirectoryFormat()

        # Call function. Patching will make sure nothing is
        # actually ran
        diamond_db = build_diamond_db(sequences)

        # Paths to inputs and outputs
        path_in = os.path.join(str(sequences), "protein-sequences.fasta")
        path_out = os.path.join(str(diamond_db), "ref_db.dmnd")

        # Check that command was called in the expected way
        cmd = [
            "diamond makedb "
            f"--in {path_in} "
            f"--db {path_out}",
            '--file-buffer-size', '67108864'
        ]

        # Check that commands is ran as expected
        subp_run.assert_called_once_with(cmd, check=True, shell=True)


class TestFetchDB(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    @patch("subprocess.run")
    def test_fetch_eggnog_db(self, subp_run):
        # Call function. Patching will make sure nothing is
        # actually ran
        eggnog_db = fetch_eggnog_db()

        # Check that command was called in the expected way
        cmd = [
            "download_eggnog_data.py", "-y", "-D",
            "--data_dir", str(eggnog_db)
        ]
        subp_run.assert_called_once_with(cmd, check=True)
