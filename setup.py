# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer

setup(
    name='q2-moshpit',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    author="Michal Ziemski, Keegan Evans",
    author_email="ziemski.michal@gmail.com",
    description="QIIME 2 plugin for metagenome analysis.",
    url="https://github.com/bokulich-lab/q2-moshpit",
    entry_points={
        'qiime2.plugins':
        ['q2-moshpit=q2_moshpit.plugin_setup:plugin']
    },
    package_data={
        'q2_moshpit': ['citations.bib', 'tests/data/*'],
        'q2_moshpit.usage_examples': ['tests/data/*'],
        'q2_moshpit.metabat2.tests': [
            'data/*', 'data/bins/*/*', 'data/contigs/*',
            'data/depth/*', 'data/maps/*', 'data/bins-small/*/*'
        ],
        'q2_moshpit.checkm.tests': [
            'data/*', 'data/bins/*', 'data/bins/*/*',
            'data/checkm_reports/*/*/*', 'data/plots/*/*/*'
        ],
        'q2_moshpit.eggnog': [
            'tests/data/*',
            'tests/data/contig-sequences-1/*',
            'tests/data/random-db-1/*',
            'tests/data/good_hits/*',
            'tests/data/bad_hits/*',
            'tests/data/eggnog_db/*',
            'tests/data/expected/*',
            'citations.bib',
        ],
        'q2_moshpit.kraken2.tests': [
            'data/*', 'data/mags/*', 'data/mags/*/*',
            'data/single-end/*', 'data/paired-end/*',
            'data/db/*', 'data/reports-mags/*',
            'data/reports-mags/*/*', 'data/outputs-mags/*',
            'data/outputs-mags/*/*', 'data/seqs/*',
            'data/bracken-db/*', 'data/bracken-report/*',
            'data/bracken-report/*/*', 'data/kraken2-report/*/*',
            'data/kraken2-reports-select/*',
            'data/kraken2-reports-select/*/*/*',
            'data/kraken-to-ncbi-tree/*',
            'data/kraken-to-ncbi-tree/*/*/*',
            'data/mags-derep/*'
        ]
    },
    zip_safe=False,
)
