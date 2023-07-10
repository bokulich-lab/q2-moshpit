# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types.per_sample_sequences import (
    SequencesWithQuality, PairedEndSequencesWithQuality
)

import q2_moshpit.kraken2.classification
from q2_types_genomics.feature_map import FeatureMap, MAGtoContigs
from q2_types_genomics.kraken2 import (
    Kraken2Reports, Kraken2Outputs, Kraken2DB
)

from qiime2.plugin import (Plugin, Citations)

import importlib
import q2_moshpit

from q2_types.sample_data import SampleData
from q2_types.feature_table import FeatureTable, Frequency
from q2_types.feature_data import FeatureData, Sequence

from q2_types_genomics.reference_db import ReferenceDB, Diamond, Eggnog
from q2_types_genomics.feature_data import NOG
from q2_types_genomics.genome_data import BLAST6
from q2_types_genomics.kraken2._type import BrackenDB
from q2_types_genomics.per_sample_data import MAGs, Contigs
from q2_types_genomics.per_sample_data._type import AlignmentMap

from qiime2.core.type import Bool, Range, Int, Str, Float, List, Choices


citations = Citations.load('citations.bib', package='q2_moshpit')

plugin = Plugin(
    name='moshpit',
    version=q2_moshpit.__version__,
    website="https://github.com/bokulich-lab/q2-moshpit",
    package='q2_moshpit',
    description=(
        'MOdular SHotgun metagenome Pipelines with Integrated '
        'provenance Tracking: QIIME 2 plugin gor metagenome analysis with'
        'tools for genome binning and functional annotation.'),
    short_description='QIIME 2 plugin for metagenome analysis.',
)

importlib.import_module('q2_moshpit.eggnog')
importlib.import_module('q2_moshpit.metabat2')

plugin.methods.register_function(
    function=q2_moshpit.metabat2.bin_contigs_metabat,
    inputs={
        'contigs': SampleData[Contigs],
        'alignment_maps': SampleData[AlignmentMap]
    },
    parameters={
        'min_contig': Int % Range(1500, None),
        'max_p': Int % Range(1, 100),
        'min_s': Int % Range(1, 100),
        'max_edges': Int % Range(1, None),
        'p_tnf': Int % Range(0, 100),
        'no_add': Bool,
        'min_cv': Int % Range(1, None),
        'min_cv_sum': Int % Range(1, None),
        'min_cls_size': Int % Range(1, None),
        'num_threads': Int % Range(0, None),
        'seed': Int % Range(0, None),
        'debug': Bool,
        'verbose': Bool
    },
    outputs=[
        ('mags', SampleData[MAGs]),
        ('contig_map', FeatureMap[MAGtoContigs])
    ],
    input_descriptions={
        'contigs': 'Placeholder.',
        'alignment_maps': 'Placeholder.'
    },
    parameter_descriptions={
        'min_contig': 'Minimum size of a contig for binning.',
        'max_p': 'Percentage of "good" contigs considered for binning '
                 'decided by connection among contigs. The greater, the '
                 'more sensitive.',
        'min_s': 'Minimum score of a edge for binning. The greater, the '
                 'more specific.',
        'max_edges': 'Maximum number of edges per node. The greater, the '
                     'more sensitive.',
        'p_tnf': 'TNF probability cutoff for building TNF graph. Use it to '
                 'skip the preparation step. (0: auto)',
        'no_add': 'Turning off additional binning for lost or small contigs.',
        'min_cv': 'Minimum mean coverage of a contig in each library '
                  'for binning.',
        'min_cv_sum': 'Minimum total effective mean coverage of a contig '
                      '(sum of depth over minCV) for binning.',
        'min_cls_size': 'Minimum size of a bin as the output.',
        'num_threads': 'Number of threads to use (0: use all cores).',
        'seed': 'For exact reproducibility. (0: use random seed)',
        'debug': 'Debug output.',
        'verbose': 'Verbose output.'
    },
    output_descriptions={
        'mags': 'The resulting MAGs.',
        'contig_map': 'MAG-to-contigs map.'
    },
    name='Bin contigs into MAGs using MetaBAT 2.',
    description='This method uses MetaBAT 2 to bin provided contigs '
                'into MAGs.',
    citations=[citations["kang2019"]]
)

plugin.methods.register_function(
    function=q2_moshpit.kraken2.classification.classify_kraken,
    inputs={
        "seqs": SampleData[
            SequencesWithQuality | PairedEndSequencesWithQuality | MAGs
            ],
        "db": Kraken2DB
    },
    parameters={
        'threads': Int % Range(1, None),
        'confidence': Float % Range(0, 1, inclusive_end=True),
        'minimum_base_quality': Int % Range(0, None),
        'memory_mapping': Bool,
        'minimum_hit_groups': Int % Range(1, None),
        'quick': Bool
    },
    outputs=[
        ('reports', SampleData[Kraken2Reports]),
        ('outputs', SampleData[Kraken2Outputs])
    ],
    input_descriptions={
        "seqs": "Sequences to be classified. Both, single-/paired-end reads"
                "and assembled MAGs, can be provided.",
        "db": "Kraken 2 database.",
    },
    parameter_descriptions={
        'threads': 'Number of threads',
        'confidence': 'Confidence score threshold.',
        'minimum_base_quality': 'Minimum base quality used in classification.'
                                ' Only applies when reads are used as input.',
        'memory_mapping': 'Avoids loading the database into RAM.',
        'minimum_hit_groups': 'Minimum number of hit groups (overlapping '
                              'k-mers sharing the same minimizer).',
        'quick': 'Quick operation (use first hit or hits).'
    },
    output_descriptions={
        'reports': 'Reports produced by Kraken2.',
        'outputs': 'Outputs produced by Kraken2.'
    },
    name='Perform taxonomic classification of bins or reads using Kraken 2.',
    description='This method uses Kraken 2 to classify provided NGS reads '
                'or MAGs into taxonomic groups.',
    citations=[citations["wood2019"]]
)

plugin.methods.register_function(
    function=q2_moshpit.kraken2.build_kraken_db,
    inputs={
        "seqs": List[FeatureData[Sequence]]
    },
    parameters={
        'collection': Str % Choices(
            ['viral', 'minusb', 'standard', 'standard8',
             'standard16', 'pluspf', 'pluspf8', 'pluspf16',
             'pluspfp', 'pluspfp8', 'pluspfp16', 'eupathdb'],
        ),
        'threads': Int % Range(1, None),
        'kmer_len': Int % Range(1, None),
        'minimizer_len': Int % Range(1, None),
        'minimizer_spaces': Int % Range(1, None),
        'no_masking': Bool,
        'max_db_size': Int % Range(0, None),
        'use_ftp': Bool,
        'load_factor': Float % Range(0, 1),
        'fast_build': Bool,
        'read_len': List[Int % Range(1, None)],
    },
    outputs=[
        ('kraken2_database', Kraken2DB),
        ('bracken_database', BrackenDB),
    ],
    input_descriptions={
        "seqs": "Sequences to be added to the Kraken 2 database."
    },
    parameter_descriptions={
        'collection': 'Name of the database collection to be fetched. '
                      'Please check https://benlangmead.github.io/aws-'
                      'indexes/k2 for the description of the available '
                      'options.',
        'threads': 'Number of threads. Only applicable when building a '
                   'custom database.',
        'kmer_len': 'K-mer length in bp/aa.',
        'minimizer_len': 'Minimizer length in bp/aa.',
        'minimizer_spaces': 'Number of characters in minimizer that are '
                            'ignored in comparisons.',
        'no_masking': 'Avoid masking low-complexity sequences prior to '
                      'building; masking requires dustmasker or segmasker '
                      'to be installed in PATH',
        'max_db_size': 'Maximum number of bytes for Kraken 2 hash table; '
                       'if the estimator determines more would normally be '
                       'needed, the reference library will be downsampled '
                       'to fit.',
        'use_ftp': 'Use FTP for downloading instead of RSYNC.',
        'load_factor': 'Proportion of the hash table to be populated.',
        'fast_build': 'Do not require database to be deterministically '
                      'built when using multiple threads. This is faster, '
                      'but does introduce variability in minimizer/LCA pairs.',
        'read_len': 'Ideal read lengths to be used while building the Bracken '
                    'database.'
    },
    output_descriptions={
        'kraken2_database': 'Kraken2 database.',
        'bracken_database': 'Bracken database.'
    },
    name='Build Kraken 2 database.',
    description='This method builds a Kraken 2/Bracken databases from '
                'provided DNA sequences or simply fetches pre-built '
                'versions from an online resource.',
    citations=[citations["wood2019"], citations["lu2017"]]
)

plugin.methods.register_function(
        function=q2_moshpit.eggnog.eggnog_diamond_search,
        inputs={'input_sequences': SampleData[Contigs],
                'diamond_db': ReferenceDB[Diamond],
                },
        parameters={
                'num_cpus': Int,
                'db_in_memory': Bool,
                },
        input_descriptions={
            'input_sequences': 'Sequence data of the contigs we want to '
                               'search for hits using the Diamond Database',
            'diamond_db': 'The filepath to an artifact containing the'
                          'Diamond database',
            },
        parameter_descriptions={
            'num_cpus': 'Number of CPUs to utilize. \'0\' will '
                        'use all available.',
            'db_in_memory': 'Read database into memory. The '
                            'database can be very large, so this '
                            'option should only be used on clusters or other '
                            'machines with enough memory.',
            },
        outputs=[('eggnog_hits', SampleData[BLAST6]),
                 ('table', FeatureTable[Frequency])
                 ],
        name='Run eggNOG search using diamond aligner',
        description="This method performs the steps by which we find our "
                    "possible target sequences to annotate using the diamond "
                    "search functionality from the eggnog `emapper.py` script",
        )

plugin.methods.register_function(
        function=q2_moshpit.eggnog.eggnog_annotate,
        inputs={
            'eggnog_hits': SampleData[BLAST6],
            'eggnog_db': ReferenceDB[Eggnog],
            },
        parameters={
            'db_in_memory': Bool,
            },
        parameter_descriptions={
            'db_in_memory': 'Read eggnog database into memory. The '
                            'eggnog database is very large(>44GB), so this '
                            'option should only be used on clusters or other '
                            'machines with enough memory.',
            },
        outputs=[('ortholog_annotations', FeatureData[NOG])],
        name='Annotate orthologs against eggNOG database',
        description="Apply eggnog mapper to annotate seed orthologs.",
        )
