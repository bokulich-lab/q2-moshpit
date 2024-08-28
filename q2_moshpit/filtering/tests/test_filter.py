# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import filecmp
import os
import shutil
import unittest
from unittest.mock import Mock, patch, ANY, call, MagicMock

import pandas as pd
import qiime2

from q2_moshpit.filtering.filter_pangenome import (
    _fetch_and_extract_grch38, _extract_fasta_from_gfa,
    _fetch_and_extract_pangenome, filter_reads_pangenome, EBI_SERVER_URL
)
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.busco.types import BUSCOResultsFormat
from q2_moshpit.filtering.filter_mags import (
    _filter_ids, _filter_manifest, _mags_to_df,
    filter_derep_mags, filter_mags
)
from q2_types.feature_data_mag import MAGSequencesDirFmt
from q2_types.per_sample_sequences import MultiMAGSequencesDirFmt


class TestMAGFiltering(TestPluginBase):
    package = 'q2_moshpit.filtering.tests'

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.metadata_df = pd.DataFrame({
            'sample_id': ['id1', 'id1', 'id2', 'id3', 'id3', 'id3'],
            'complete': [20.0, 98.0, 68.5, 100.0, 98.0, 21.3]
        }, index=pd.Index(
            ['mag1', 'mag2', 'mag3', 'mag4', 'mag5', 'mag6'],
            name='id'
        ))

    def setUp(self):
        super().setUp()
        self.manifest_df = pd.read_csv(self.get_data_path("MANIFEST"))
        self.manifest_df.set_index(['sample-id', 'mag-id'], inplace=True)
        self.mag_data_dir = self.get_data_path('mags')
        self.mag_derep_data_dir = self.get_data_path('mags/sample2')
        self.mag_df = pd.DataFrame({
            'sample_id': ['sample1', 'sample2', 'sample2'],
            'mag_id': [
                '24dee6fe-9b84-45bb-8145-de7b092533a1',
                'd65a71fa-4279-4588-b937-0747ed5d604d',
                'db03f8b6-28e1-48c5-a47c-9c65f38f7357'
            ],
            'mag_fp': [
                f'{self.mag_data_dir}/sample1/'
                f'24dee6fe-9b84-45bb-8145-de7b092533a1.fasta',
                f'{self.mag_data_dir}/sample2/'
                f'd65a71fa-4279-4588-b937-0747ed5d604d.fasta',
                f'{self.mag_data_dir}/sample2/'
                f'db03f8b6-28e1-48c5-a47c-9c65f38f7357.fasta'
            ]
        })
        self.metadata = qiime2.Metadata(self.metadata_df)

    def test_filter_ids_include(self):
        ids = {'mag1', 'mag2', 'mag3', 'mag4', 'mag5', 'mag6'}
        obs = _filter_ids(
            ids, self.metadata, 'complete>40', exclude_ids=False
        )
        exp = {'mag2', 'mag3', 'mag4', 'mag5'}
        self.assertSetEqual(obs, exp)

    def test_filter_ids_include_none(self):
        ids = {'mag1', 'mag2', 'mag3', 'mag4', 'mag5', 'mag6'}
        obs = _filter_ids(
            ids, self.metadata, 'complete<10', exclude_ids=False
        )
        exp = ids
        self.assertSetEqual(obs, exp)

    def test_filter_ids_exclude(self):
        ids = {'mag1', 'mag2', 'mag3', 'mag4', 'mag5', 'mag6'}
        obs = _filter_ids(
            ids, self.metadata, 'complete>40', exclude_ids=True
        )
        exp = {'mag1', 'mag6'}
        self.assertSetEqual(obs, exp)

    def test_filter_manifest_mags(self):
        ids = {'mag1', 'mag2', 'mag5'}
        obs = _filter_manifest(self.manifest_df, ids, on='mag')
        exp = self.manifest_df[
            self.manifest_df.index.get_level_values('mag-id').isin(ids)
        ]
        pd.testing.assert_frame_equal(obs, exp)

    def test_filter_manifest_samples(self):
        ids = {'id1', }
        obs = _filter_manifest(self.manifest_df, ids, on='sample')
        exp = self.manifest_df[
            self.manifest_df.index.get_level_values('sample-id').isin(ids)
        ]
        pd.testing.assert_frame_equal(obs, exp)

    def test_filter_manifest_error(self):
        ids = {'id1', }
        with self.assertRaisesRegex(ValueError, 'parameter: unicorn'):
            _filter_manifest(self.manifest_df, ids, on='unicorn')

    def test_mags_to_df_on_sample(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode='r')
        obs = _mags_to_df(mags, on='sample')
        exp = self.mag_df
        exp.set_index('sample_id', inplace=True)
        pd.testing.assert_frame_equal(obs, exp)

    def test_mags_to_df_on_mag(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode='r')
        obs = _mags_to_df(mags, on='mag')
        exp = self.mag_df
        exp.set_index('mag_id', inplace=True)
        pd.testing.assert_frame_equal(obs, exp)

    def test_filter_derep_mags(self):
        mags = MAGSequencesDirFmt(self.mag_derep_data_dir, mode='r')
        metadata = BUSCOResultsFormat(
            self.get_data_path('metadata-derep.tsv'), mode='r'
        ).view(qiime2.Metadata)

        obs = filter_derep_mags(mags, metadata, where='complete>10')
        obs_features = obs.feature_dict()
        exp_features = ['db03f8b6-28e1-48c5-a47c-9c65f38f7357']
        self.assertListEqual(list(obs_features.keys()), exp_features)

    def test_filter_mags_features(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode='r')
        metadata = BUSCOResultsFormat(
            self.get_data_path('metadata-derep.tsv'), mode='r'
        ).view(qiime2.Metadata)

        obs = filter_mags(mags, metadata, where='length<1000000', on='mag')
        obs_samples = obs.sample_dict()
        exp_samples = ['sample2']
        self.assertListEqual(list(obs_samples.keys()), exp_samples)

        obs_features = [
            mag_id for feature_dict in obs_samples.values()
            for mag_id in feature_dict.keys()
        ]
        exp_features = [
            'd65a71fa-4279-4588-b937-0747ed5d604d',
            'db03f8b6-28e1-48c5-a47c-9c65f38f7357'
        ]
        self.assertListEqual(obs_features, exp_features)

    def test_filter_mags_samples(self):
        mags = MultiMAGSequencesDirFmt(self.mag_data_dir, mode='r')
        metadata = qiime2.Metadata(
            pd.read_csv(self.get_data_path('metadata-sample.tsv'),
                        sep='\t', index_col=0)
        )

        obs = filter_mags(mags, metadata, where='metric<5', on='sample')
        obs_samples = obs.sample_dict()
        exp_samples = ['sample2']
        self.assertListEqual(list(obs_samples.keys()), exp_samples)

        obs_feature_count = len(obs.sample_dict()['sample2'])
        exp_feature_count = 2
        self.assertEqual(obs_feature_count, exp_feature_count)

    @patch('shutil.move')
    def test_fetch_and_extract_grch38(self, p1):
        fake_results = Mock()
        fake_callable = Mock(return_value=fake_results)

        _fetch_and_extract_grch38(fake_callable, "/some/where")

        fake_callable.assert_called_once_with(
            taxon='Homo sapiens',
            only_reference=True,
            assembly_levels=['chromosome'],
            assembly_source='refseq',
            only_genomic=False
        )
        fake_results.genome_assemblies.export_data.assert_called_once_with(
            "/some/where"
        )
        p1.assert_called_once_with(
            "/some/where/dna-sequences.fasta", "/some/where/grch38.fasta"
        )

    @patch("subprocess.run")
    @patch("os.remove")
    def test_extract_from_gfa(self, p1, p2):
        fasta_fp = os.path.join(self.temp_dir.name, "some_fasta.fa")
        _extract_fasta_from_gfa("/some/gfa", fasta_fp)

        p2.assert_called_once_with(
            ["gfatools", "gfa2fa", "/some/gfa"], stdout=ANY
        )
        p1.assert_called_once_with("/some/gfa")

    @patch("subprocess.run", side_effect=OSError)
    def test_extract_from_gfa_error(self, p1):
        fasta_fp = os.path.join(self.temp_dir.name, "some_fasta.fa")
        with self.assertRaisesRegex(
                Exception, "Failed to extract"
        ):
            _extract_fasta_from_gfa("/some/gfa", fasta_fp)

    @patch("q2_moshpit.filtering.filter_pangenome.run_command")
    def test_fetch_and_extract_pangenome(self, p1):
        uri = "http://hello.org/file123.gz"
        _fetch_and_extract_pangenome(uri, "/some/where")

        p1.assert_has_calls([
            call(["wget", uri, "-q", "-O", "/some/where/file123.gz"]),
            call(["gunzip", "/some/where/file123.gz"]),
        ])

    @patch(
        "q2_moshpit.filtering.filter_pangenome.run_command",
        side_effect=OSError
    )
    def test_fetch_and_extract_pangenome_error(self, p1):
        with self.assertRaisesRegex(Exception, "Unable to connect"):
            _fetch_and_extract_pangenome("http://hello.org", "/some/where")

    @patch(
        'q2_moshpit.filtering.filter_pangenome._fetch_and_extract_pangenome'
    )
    @patch('q2_moshpit.filtering.filter_pangenome._fetch_and_extract_grch38')
    @patch('q2_moshpit.filtering.filter_pangenome._extract_fasta_from_gfa')
    def test_filter_reads_pangenome(
            self, mock_extract_fasta, mock_fetch_grch38, mock_fetch_pangenome
    ):
        ctx = MagicMock()
        ctx.get_action.return_value = MagicMock()

        mock_build_index_result = MagicMock()
        ctx.get_action(
            "quality_control", "bowtie2_build"
        ).return_value = (mock_build_index_result,)

        mock_filtered_reads_result = MagicMock()
        ctx.get_action(
            "quality_control", "filter_reads"
        ).return_value = (mock_filtered_reads_result,)
        ctx.make_artifact.return_value = MagicMock()

        reads = MagicMock()
        index = None
        n_threads = 1

        open(os.path.join(self.temp_dir.name, "pangenome.gfa"), 'w').close()
        shutil.copy(
            self.get_data_path("pangenome/grch38.fasta"),
            os.path.join(self.temp_dir.name, "grch38.fasta")
        )
        shutil.copy(
            self.get_data_path("pangenome/pangenome.fasta"),
            os.path.join(self.temp_dir.name, "pangenome.fasta")
        )
        with patch(
                'tempfile.TemporaryDirectory',
                return_value=self.temp_dir
        ) as mock_tempdir:
            filtered_reads, generated_index = filter_reads_pangenome(
                ctx=ctx,
                reads=reads,
                index=index,
                n_threads=n_threads
            )

            # Assertions
            ctx.get_action.assert_any_call("rescript", "get_ncbi_genomes")
            ctx.get_action.assert_any_call("quality_control", "bowtie2_build")
            ctx.get_action.assert_any_call("quality_control", "filter_reads")

            mock_fetch_pangenome.assert_called_once_with(
                EBI_SERVER_URL, self.temp_dir.name
            )
            mock_fetch_grch38.assert_called_once_with(
                ctx.get_action("rescript", "get_ncbi_genomes"),
                self.temp_dir.name
            )
            mock_extract_fasta.assert_called_once_with(
                os.path.join(self.temp_dir.name, "pangenome.gfa"),
                os.path.join(self.temp_dir.name, "pangenome.fasta")
            )
            assert filecmp.cmp(
                self.get_data_path("pangenome/combined.fasta"),
                os.path.join(self.temp_dir.name, "combined.fasta")
            )

            ctx.make_artifact.assert_called_once(
                "FeatureData[Sequence]",
                os.path.join(self.temp_dir.name, "combined.fasta")
            )
            ctx.get_action(
                "quality_control", "bowtie2_build"
            ).assert_called_once()
            ctx.get_action(
                "quality_control", "filter_reads"
            ).assert_called_once_with(
                demultiplexed_sequences=reads,
                database=index,
                exclude_seqs=True,
                n_threads=n_threads
            )
            self.assertIsNotNone(generated_index)


if __name__ == "__main__":
    unittest.main()
