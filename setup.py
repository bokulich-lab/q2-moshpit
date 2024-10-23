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
    author='Michal Ziemski, Keegan Evans',
    author_email='ziemski.michal@gmail.com',
    description='QIIME 2 plugin for metagenome analysis.',
    url='https://github.com/bokulich-lab/q2-moshpit',
    entry_points={
        'qiime2.plugins':
        ['q2-moshpit=q2_moshpit.plugin_setup:plugin']
    },
    package_data={
        '.': ['.coveragerc',],
        'q2_moshpit': [
            'citations.bib',
            'tests/data/*',
            'tests/data/mags-derep/*',
            'tests/data/md5/*',
            'assets/busco/*',
            'assets/busco/feature_data/*',
            'assets/busco/sample_data/*',
            'assets/busco/js/*',
            'assets/busco/css/*',
        ],
        'q2_moshpit.usage_examples': [
            'tests/data/*'
        ],
        'q2_moshpit.metabat2.tests': [
            'data/*',
            'data/bins/*/*',
            'data/contigs/*',
            'data/depth/*',
            'data/maps/*',
            'data/bins-small/*/*',
            'data/bins-no-uuid/*/*'
        ],
        'q2_moshpit.busco.tests': [
            'data/*',
            'data/summaries/*',
            'data/mags/*',
            'data/mags/sample1/*',
            'data/mags/sample2/*',
            'data/busco_db/busco_downloads/*',
            'data/busco_db/busco_downloads/placement_files/*',
            'data/busco_db/busco_downloads/lineages/lineage_1/*',
            'data/busco_db/busco_downloads/lineages/lineage_1/hmms/*',
            'data/busco_db/busco_downloads/lineages/lineage_1/info/*',
            'data/busco_db/busco_downloads/lineages/lineage_1/prfl/*'
            'data/busco_results/*',
            'data/busco_results/sample1/*',
            'data/busco_results/sample2/*',
            'data/busco_results/collated/*',

        ],
        'q2_moshpit.busco.types.tests': [
            'data/*',
            'data/*/*',
            'data/*/*/*',
            'data/*/*/*/*',
            'data/*/*/*/*/*',
            'data/*/*/*/*/*/*',
        ],
        'q2_moshpit.eggnog.tests': [
            'data/*',
            'data/annotations/*',
            'data/build_eggnog_diamond_db/*',
            'data/contig-sequences-1/*',
            'data/eggnog_db/*',
            'data/expected/*',
            'data/fastas/*',
            'data/good_hits/*',
            'data/hits/*',
            'data/hmmer/fastas/1/*',
            'data/hmmer/hmms/1/*',
            'data/idmap/*',
            'data/invalid_idmaps/*',
            'data/mag-sequences/*',
            'data/mag-sequences-per-sample/*',
            'data/mag-sequences-per-sample/sample1/*',
            'data/mag-sequences-per-sample/sample2/*',
            'data/md5/*',
            'data/pressed_hmm/*',
            'data/random-db-1/*',
            'data/valid_idmaps/*',
            'citations.bib',
        ],
        'q2_moshpit.kraken2.tests': [
            'data/*',
            'data/mags-derep/*',
            'data/mags/*/*',
            'data/single-end/*',
            'data/paired-end/*',
            'data/db/*',
            'data/reports-mags/*',
            'data/reports-mags/*/*',
            'data/reports-mags-unclassified-missing-frac/*',
            'data/reports-mags-unclassified-wrong-frac/*',
            'data/outputs-mags/*',
            'data/outputs-mags/*/*',
            'data/seqs/*',
            'data/bracken-db/*',
            'data/bracken-report/*',
            'data/bracken-report/*/*',
            'data/bracken-report-with-unclassified/*',
            'data/bracken-report-with-unclassified/*/*',
            'data/kraken2-reports-select/*',
            'data/kraken2-reports-select/*/*',
            'data/kraken2-to-ncbi-tree/*',
            'data/kraken2-to-ncbi-tree/*/*/*',
            'data/simulated-sequences/**/*',
        ],
        'q2_moshpit.dereplication.tests': [
            'data/*',
            'data/mags/*',
            'data/mags/*/*',
            'data/mags-unique/*',
        ],
        'q2_moshpit.kaiju.tests': [
            'data/*',
            'data/*/*'
        ],
        'q2_moshpit.prodigal.tests': [
            'data/*',
            'data/*/*',
        ],
        'q2_moshpit.partition.tests': [
            'data/*',
            'data/*/*',
            'data/*/*/*',
            'data/*/*/*/*',
            'data/*/*/*/*/*',
        ],
        'q2_moshpit.filtering.tests': [
            'data/*',
            'data/mags/*',
            'data/mags/*/*',
            'data/pangenome/*'
        ],
        'q2_moshpit.abundance.tests': [
            'data/*',
            'data/mag-index/*',
            'data/mag-sequences/*',
            'data/reads/*',
            'data/reads-mapped/*'
        ]
    },
    zip_safe=False,
)
