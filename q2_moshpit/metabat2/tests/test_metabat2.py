# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob

import shutil

import os
import tempfile
import unittest
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt, BAMDirFmt
from q2_types_genomics.per_sample_data._format import MultiFASTADirectoryFormat
from unittest.mock import patch, ANY, call

from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.metabat2.metabat2 import (_assert_samples,
                                          _get_sample_name_from_path,
                                          _sort_bams, _estimate_depth,
                                          _run_metabat2, _process_sample,
                                          _bin_contigs_metabat, _generate_contig_map)


class TestMetabat2(TestPluginBase):
    package = 'q2_moshpit.metabat2.tests'

    def test_get_sample_name_from_path(self):
        obs = _get_sample_name_from_path('/a/b/sampleX.fasta',
                                         '.fasta')
        exp = 'sampleX'
        self.assertEqual(exp, obs)

    def test_get_sample_name_from_path_underscores(self):
        obs = _get_sample_name_from_path('/a/b/sampleX_something.fasta',
                                         '_something.fasta')
        exp = 'sampleX'
        self.assertEqual(exp, obs)

    def test_assert_samples_ok(self):
        contigs = ['/a/b/s1_contigs.fa', '/a/b/s3_contigs.fa',
                   '/a/b/s2_contigs.fa']
        maps = ['/a/b/s3_alignment.bam', '/a/b/s2_alignment.bam',
                '/a/b/s1_alignment.bam']

        obs_samples = _assert_samples(contigs, maps)
        exp_samples = {
            's1': {'contigs': '/a/b/s1_contigs.fa',
                   'map': '/a/b/s1_alignment.bam'},
            's2': {'contigs': '/a/b/s2_contigs.fa',
                   'map': '/a/b/s2_alignment.bam'},
            's3': {'contigs': '/a/b/s3_contigs.fa',
                   'map': '/a/b/s3_alignment.bam'}
        }
        self.assertDictEqual(exp_samples, obs_samples)

    def test_assert_samples_uneven(self):
        contigs = ['/a/b/s1_contigs.fa', '/a/b/s3_contigs.fa']
        maps = ['/a/b/s3_alignment.bam', '/a/b/s2_alignment.bam',
                '/a/b/s1_alignment.bam']

        with self.assertRaisesRegex(
                Exception,
                'Contigs and alignment maps should belong to the same sample'
                ' set. You provided contigs for samples: s1,s3 but maps for'
                ' samples: s1,s2,s3. Please check your inputs and try again.'
        ):
            _assert_samples(contigs, maps)

    def test_assert_samples_non_matching(self):
        contigs = ['/a/b/s1_contigs.fa', '/a/b/s4_contigs.fa',
                   '/a/b/s2_contigs.fa']
        maps = ['/a/b/s3_alignment.bam', '/a/b/s2_alignment.bam',
                '/a/b/s1_alignment.bam']

        with self.assertRaisesRegex(
                Exception,
                'contigs for samples: s1,s2,s4 but maps for samples: s1,s2,s3'
        ):
            _assert_samples(contigs, maps)

    @patch('subprocess.run')
    def test_sort_bams_ok(self, p1):
        fake_props = {'map': '/some/where/map.bam', 'fake_key': 'abc'}

        obs_props = _sort_bams('samp1', fake_props, '/new/location')
        exp_props = {
            'map': '/new/location/samp1_alignment_sorted.bam',
            'fake_key': 'abc'
        }

        self.assertDictEqual(exp_props, obs_props)
        p1.assert_called_once_with(
            ['samtools', 'sort', fake_props['map'], '-o',
             '/new/location/samp1_alignment_sorted.bam'], check=True
        )

    @patch('subprocess.run')
    def test_estimate_depth_ok(self, p1):
        fake_props = {'map': '/some/where/map.bam', 'fake_key': 'abc'}

        obs_fp = _estimate_depth('samp1', fake_props, '/new/location')
        exp_fp = '/new/location/samp1_depth.txt'

        self.assertEqual(exp_fp, obs_fp)
        p1.assert_called_once_with(
            ['jgi_summarize_bam_contig_depths', '--outputDepth',
             exp_fp, fake_props['map']], check=True
        )

    @patch('subprocess.run')
    def test_run_metabat2_ok(self, p1):
        fake_props = {'map': '/some/where/map.bam', 'contigs': '/a/b/co.fa'}
        fake_args = ['--verbose', '--pTNF', '85', '--minContig', '2500']

        with tempfile.TemporaryDirectory() as fake_loc:
            obs_fp = _run_metabat2('samp1', fake_props, fake_loc,
                                   '/some/depth/file.txt', fake_args)
            exp_fp = os.path.join(fake_loc, 'samp1')

            self.assertEqual(exp_fp, obs_fp)

            exp_cmd = [
                'metabat2', '-i', fake_props['contigs'],
                '-a', '/some/depth/file.txt', '-o',
                os.path.join(fake_loc, 'samp1', 'bin')
            ]
            exp_cmd.extend(fake_args)
            p1.assert_called_once_with(exp_cmd, check=True)

    @patch('q2_moshpit.metabat2.uuid4')
    @patch('q2_moshpit.metabat2._sort_bams')
    @patch('q2_moshpit.metabat2._estimate_depth')
    @patch('q2_moshpit.metabat2._run_metabat2')
    def test_process_sample(self, p1, p2, p3, p4):
        fake_props = {
            'map': 'some/where/samp1_alignment.bam',
            'contigs': 'some/where/samp1_contigs.fasta'
        }
        fake_props_mod = {
            'map': 'some/where/samp1_alignment_sorted.bam',
            'contigs': 'some/where/samp1_contigs.fasta'
        }
        fake_args = ['--verbose', '--minContig', '1500', '--minClsSize',
                     '10000']

        p3.return_value = fake_props_mod
        p2.return_value = 'some/where/samp1_depth.txt'
        p4.side_effect = [
            '522775d4-b1c6-4ee3-8b47-cd990f17eb8b',
            '684db670-6304-4f33-a0ea-7f570532e178',
            '37356c23-b8db-4bbe-b4c9-d35e1cef615b',
            '51c19113-31f0-4e4c-bbb3-b9df26b949f3'
        ]

        with tempfile.TemporaryDirectory() as fake_loc:
            p1.return_value = os.path.join(fake_loc, 'bins', 'samp1')

            # copy two expected bins to the new location
            samp1_bins_fp = self.get_data_path('bins/samp1')
            shutil.copytree(
                samp1_bins_fp, os.path.join(fake_loc, 'bins', 'samp1'),
                dirs_exist_ok=True
            )

            _process_sample('samp1', fake_props, fake_args, fake_loc)

            # find the newly formed bins
            obs_bins = set([
                x.split('/')[-1] for x in
                glob.glob(os.path.join(fake_loc, 'samp1', '*.fasta'))
            ])
            exp_bins = {
                '522775d4-b1c6-4ee3-8b47-cd990f17eb8b.fasta',
                '684db670-6304-4f33-a0ea-7f570532e178.fasta'
            }
            self.assertSetEqual(exp_bins, obs_bins)

            p3.assert_called_once_with('samp1', fake_props, ANY)
            p2.assert_called_once_with('samp1', fake_props_mod, ANY)
            p1.assert_called_once_with(
                'samp1', fake_props_mod, ANY,
                'some/where/samp1_depth.txt', fake_args
            )

    @patch('q2_moshpit.metabat2.MultiFASTADirectoryFormat')
    @patch('q2_moshpit.metabat2._process_sample')
    def test_bin_contigs_metabat(self, p1, p2):
        input_contigs = self.get_data_path('contigs')
        input_maps = self.get_data_path('maps')
        contigs = ContigSequencesDirFmt(input_contigs, mode='r')
        maps = BAMDirFmt(input_maps, mode='r')

        args = ['--verbose', '--minContig', '1500', '--minClsSize', '10000']

        mock_contigs = MultiFASTADirectoryFormat(
            self.get_data_path('bins'), 'r'
        )
        p2.return_value = mock_contigs

        obs_bins, obs_map = _bin_contigs_metabat(contigs, maps, args)

        self.assertIsInstance(obs_bins, MultiFASTADirectoryFormat)
        p1.assert_has_calls([
            call(
                'samp1',
                {'contigs': self.get_data_path('/contigs/samp1_contigs.fa'),
                 'map': self.get_data_path('/maps/samp1_alignment.bam')},
                args, str(mock_contigs)
            ),
            call(
                'samp2',
                {'contigs': self.get_data_path('/contigs/samp2_contigs.fa'),
                 'map': self.get_data_path('/maps/samp2_alignment.bam')},
                args, str(mock_contigs)
            )
        ])

        # find the newly formed bins
        obs_bins = sorted([sorted([
            '/'.join(x.split('/')[-2:]) for x in
            glob.glob(os.path.join(str(obs_bins), f'samp{y}', '*.fa'))
        ]) for y in (1, 2)])
        exp_bins = [
            ['samp1/522775d4-b1c6-4ee3-8b47-cd990f17eb8b.fa',
             'samp1/684db670-6304-4f33-a0ea-7f570532e178.fa'],
            ['samp2/37356c23-b8db-4bbe-b4c9-d35e1cef615b.fa',
             'samp2/51c19113-31f0-4e4c-bbb3-b9df26b949f3.fa']
        ]
        self.assertListEqual(exp_bins, obs_bins)

    def test_generate_contig_map(self):
        contigs = MultiFASTADirectoryFormat(
            self.get_data_path('bins-small'), 'r'
        )
        obs = _generate_contig_map(contigs)
        exp = {
            '684db670-6304-4f33-a0ea-7f570532e178': [
                'NODE_2', 'NODE_6', 'NODE_7'
            ],
            '522775d4-b1c6-4ee3-8b47-cd990f17eb8b': [
                'NODE_8', 'NODE_11'
            ],
            '51c19113-31f0-4e4c-bbb3-b9df26b949f3': [
                'NODE_12', 'NODE_13', 'NODE_14'
            ],
            '37356c23-b8db-4bbe-b4c9-d35e1cef615b': [
                'NODE_2', 'NODE_6', 'NODE_7', 'NODE_15'
            ]
        }
        self.assertDictEqual(exp, obs)


if __name__ == '__main__':
    unittest.main()
