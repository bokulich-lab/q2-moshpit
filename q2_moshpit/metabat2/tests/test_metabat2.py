# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
import unittest
from unittest.mock import patch

from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.metabat2.metabat2 import (_assert_samples,
                                          _get_sample_name_from_path,
                                          _sort_bams, _estimate_depth,
                                          _run_metabat2)


class TestMetabat2(TestPluginBase):
    package = 'q2_moshpit.metabat2.tests'

    def test_get_sample_name_from_path(self):
        obs = _get_sample_name_from_path('/a/b/sampleX.fasta')
        exp = 'sampleX'
        self.assertEqual(exp, obs)

    def test_get_sample_name_from_path_underscores(self):
        obs = _get_sample_name_from_path('/a/b/sampleX_something.fasta')
        exp = 'sampleX'
        self.assertEqual(exp, obs)

    def test_assert_samples_ok(self):
        contigs = ['/a/b/s1_tigs.fa', '/a/b/s3_tigs.fa', '/a/b/s2_tigs.fa']
        maps = ['/a/b/s3_aln.bam', '/a/b/s2_aln.bam', '/a/b/s1_aln.bam']

        obs_samples = _assert_samples(contigs, maps)
        exp_samples = {
            's1': {'contigs': '/a/b/s1_tigs.fa', 'map': '/a/b/s1_aln.bam'},
            's2': {'contigs': '/a/b/s2_tigs.fa', 'map': '/a/b/s2_aln.bam'},
            's3': {'contigs': '/a/b/s3_tigs.fa', 'map': '/a/b/s3_aln.bam'}
        }
        self.assertDictEqual(exp_samples, obs_samples)

    def test_assert_samples_uneven(self):
        contigs = ['/a/b/s1_tigs.fa', '/a/b/s3_tigs.fa']
        maps = ['/a/b/s3_aln.bam', '/a/b/s2_aln.bam', '/a/b/s1_aln.bam']

        with self.assertRaisesRegex(
                Exception,
                'contigs for samples: s1,s3 but maps for samples: s1,s2,s3'
        ):
            _assert_samples(contigs, maps)

    def test_assert_samples_non_matching(self):
        contigs = ['/a/b/s1_tigs.fa', '/a/b/s4_tigs.fa', '/a/b/s2_tigs.fa']
        maps = ['/a/b/s3_aln.bam', '/a/b/s2_aln.bam', '/a/b/s1_aln.bam']

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
        fake_args = ['--verbose', '--pTNF', '85', '--min-contig', '2500']

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


if __name__ == '__main__':
    unittest.main()
