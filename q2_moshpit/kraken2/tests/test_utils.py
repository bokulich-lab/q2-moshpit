# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest

import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.kraken2.utils import _process_kraken2_arg, _find_lca, _find_lca_majority, _find_super_lca


class TestKraken2Utils(TestPluginBase):
    package = 'q2_moshpit.kraken2.tests'

    def setUp(self):
        super().setUp()
        m_flor = [
            ['Bacteria', 'Bacteria', 'Actinomycetota', 'Actinomycetes',
             'Mycobacteriales', 'Mycobacteriaceae', 'Mycobacterium',
             'Mycobacterium florentinum']
        ] * 5
        m_methano = [
            ['Bacteria', 'Bacteria', 'Actinomycetota', 'Actinomycetes',
             'Mycobacteriales', 'Mycobacteriaceae', 'Mycobacterium',
             'Candidatus Mycobacterium methanotrophicum']
        ] * 3
        self.taxa = pd.Series([
            *m_flor, *m_methano,
            ['Bacteria', 'Bacteria', 'Actinomycetota', 'Actinomycetes',
             'Mycobacteriales', 'Mycobacteriaceae', 'Mycobacterium',
             'Mycobacterium avium',
             'Mycobacterium avium subsp. hominissuis'],
            ['Bacteria', 'Bacteria', 'Actinomycetota', 'Actinomycetes',
             'Mycobacteriales', 'Mycobacteriaceae', 'Mycobacterium',
             'Mycobacterium heckeshornense'],
            ['Bacteria', 'Bacteria', 'Actinomycetota', 'Actinomycetes',
             'Mycobacteriales', 'Mycobacteriaceae', 'Mycobacterium',
             'Mycobacterium heckeshornense']
        ], name='Taxon')

    def test_process_kraken2_arg_bool(self):
        obs = _process_kraken2_arg('quick', True)
        exp = ['--quick']
        self.assertListEqual(obs, exp)

    def test_process_kraken2_arg_number(self):
        obs = _process_kraken2_arg('threads', 3)
        exp = ['--threads', '3']
        self.assertListEqual(obs, exp)

    def test_process_kraken2_arg_string(self):
        obs = _process_kraken2_arg('db', '/some/where/test_db')
        exp = ['--db', '/some/where/test_db']
        self.assertListEqual(obs, exp)

    def test_process_kraken2_arg_unknown(self):
        with self.assertRaisesRegex(
                NotImplementedError,
                'Parsing arguments of type "<class \'list\'>" '
                'is not supported.'
        ):
            _process_kraken2_arg('fake_param', [1, 2])

    def test_find_lca(self):
        obs = list(_find_lca(self.taxa))
        exp = [
            'Bacteria', 'Bacteria', 'Actinomycetota', 'Actinomycetes',
            'Mycobacteriales', 'Mycobacteriaceae', 'Mycobacterium'
        ]
        self.assertListEqual(obs, exp)

    def test_find_lca_majority(self):
        obs = list(_find_lca_majority(self.taxa))
        exp = [
            'Bacteria', 'Bacteria', 'Actinomycetota', 'Actinomycetes',
            'Mycobacteriales', 'Mycobacteriaceae', 'Mycobacterium',
            'Mycobacterium florentinum'
        ]
        self.assertListEqual(obs, exp)

    def test_find_super_lca(self):
        obs = list(_find_super_lca(self.taxa))
        exp = [
            'Bacteria', 'Bacteria', 'Actinomycetota', 'Actinomycetes',
            'Mycobacteriales', 'Mycobacteriaceae', 'Mycobacterium',
            'Mycobacterium florentinum'
        ]
        self.assertListEqual(obs, exp)


if __name__ == '__main__':
    unittest.main()
