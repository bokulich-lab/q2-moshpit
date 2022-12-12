# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import importlib
import q2_moshpit
from q2_types.sample_data import SampleData

from q2_types_genomics.reference_db import ReferenceDB, Diamond, Eggnog

from q2_types_genomics.ortholog import Ortholog, Seed, Annotation

from q2_types_genomics.per_sample_data import MAGs, Contigs
from q2_types_genomics.per_sample_data._type import AlignmentMap
from qiime2.plugin import Bool, Range, Int
from qiime2.plugin import Plugin, Citations
import q2_moshpit.usage_examples._examples as all_xmpls


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

importlib.import_module('q2_moshpit.diamond')
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
    outputs=[('mags', SampleData[MAGs])],
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
    output_descriptions={'mags': 'The resulting MAGs.'},
    name='Bin contigs into MAGs using MetaBAT 2.',
    description='This method uses MetaBAT 2 to bin provided contigs '
                'into MAGs.',
    citations=[]
)
# diamond_search
plugin.methods.register_function(
        function=q2_moshpit.diamond.eggnog_diamond_search,
        inputs={'input_sequences': SampleData[Contigs],
                'diamond_db': ReferenceDB[Diamond],
                },
        parameters={
                'num_cpus': Int,
                },
        input_descriptions={
            'input_sequences': ('Sequence data of the contigs we want to '
                                'search for hits using the Diamond Database'),
            'diamond_db': 'The filepath to an artifact containing the'
                          'Diamond database',
            },
        parameter_descriptions={
            'num_cpus': 'Number of CPUs to utilize. \'0\' will '
                        'use all available.',
            },
        outputs=[('seed_ortholog', Ortholog[Seed])],
        name='eggnog_diamond_search',
        description="This method performs the steps by which we find our "
                    "possible target sequences to annotate using the diamond "
                    "search functionality from the eggnog `emapper.py` script",
        examples={'eggnog_diamond_search':
                  all_xmpls.eggnog_diamond_search_example},
        )

plugin.methods.register_function(
        function=q2_moshpit.annotation.eggnog_annotate_seed_orthologs,
        inputs={'hits_table': Ortholog[Seed],
                'eggnog_db': ReferenceDB[Eggnog],
                },
        parameters={},
        outputs=[('annotation_ortholog', Ortholog[Annotation])],
        name='eggnog_annotate_seed_orthologs',
        description="Uses Eggnog Mapper to apply functional annotations from "
        "the eggnog database to previously generated \"seed orthologs\".",
        examples={'eggnog_annotate_seed_orthologs':
                  all_xmpls.eggnog_annotate_seed_orthologs_example},
        )
