# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest

from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.metabat2.metabat2 import _assert_samples


class TestMetabat2(TestPluginBase):
    package = 'q2_moshpit.metabat2.tests'

    def test_assert_samples_ok(self):
        contigs = ['/a/b/s1_tigs.fa', '/a/b/s3_tigs.fa', '/a/b/s2_tigs.fa']
        maps = ['/a/b/s3_aln.bam', '/a/b/s2_aln.bam', '/a/b/s1_aln.bam']

        obs_samples = _assert_samples(contigs, maps)
        exp_samples = {
            's1': {'contigs': '/a/b/s1_tigs.fa', 'map': '/a/b/s1_aln.bam'},
            's2': {'contigs': '/a/b/s2_tigs.fa', 'map': '/a/b/s2_aln.bam'},
            's3': {'contigs': '/a/b/s3_tigs.fa', 'map': '/a/b/s3_aln.bam'}
        }
        self.assertDictEqual(obs_samples, exp_samples)

    def test_assert_samples_uneven(self):
        contigs = ['/a/b/s1_tigs.fa', '/a/b/s3_tigs.fa']
        maps = ['/a/b/s3_aln.bam', '/a/b/s2_aln.bam', '/a/b/s1_aln.bam']

        with self.assertRaisesRegex(
            Exception, 'provided contigs for 2 samples.*maps for 3 samples.'
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


if __name__ == '__main__':
    unittest.main()
