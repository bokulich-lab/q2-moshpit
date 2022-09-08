# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.sample_data import SampleData
from q2_types.feature_data import FeatureData, Sequence

from q2_types_genomics.feature_data import (
    DiamondDB, MMseq2DB, NOG
    )
from q2_types_genomics.eggnog import ReferenceDB, Eggnog
from q2_types_genomics.per_sample_data import MAGs, Contigs
from q2_types_genomics.per_sample_data._type import AlignmentMap
from qiime2.plugin import Bool, Range, Int, Str, Choices, List
from qiime2.plugin import (Plugin, Citations, TypeMap)

import q2_moshpit

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


T_mode, T_OUT_fmt = TypeMap({
    Str % Choices("diamond"): FeatureData[DiamondDB],
    Str % Choices("mmseqs"): FeatureData[MMseq2DB],
})

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

plugin.methods.register_function(
    function=q2_moshpit.eggnog.create_reference_db,
    inputs={},
    parameters={'mode': T_mode,
                'target_taxa': List[Str],
                'name': Str,
                'simulate': Bool,
                },
    outputs=[('ref_db', T_OUT_fmt),
             ],
    name='download_diamond_db',
    description='Uses EggnogMapper\'s built in utility to download'
                'a Diamond reference database',
)

plugin.methods.register_function(
        function=q2_moshpit.eggnog.e_mapper,
        inputs={'main_db': ReferenceDB[Eggnog],
                'ancillary_db': FeatureData[DiamondDB | MMseq2DB],
                'query_sequences': FeatureData[Sequence],
                },
        parameters={'itype': Str % Choices(["metagenome", "genome", "CDS",
                                            "proteins",
                                            ]),
                    'result_name': Str,
                    'mode': Str % Choices(["diamond", "mmseqs"]),
                    },
        outputs=[('annotation_file', FeatureData[NOG]),
                 ],
        name='Create annotation mappings',
        description='Annotatates target features with taxonomic data',
)

plugin.methods.register_function(
        function=q2_moshpit.eggnog.download_references,
        inputs={},
        parameters={'simulate': Bool,
                    'target_taxa': List[Str],
                    },
        outputs=[('references', ReferenceDB[Eggnog]),
                 ],
        name='download_references',
        description='Downloads required reference databases for performing'
                    'annotations using eggnog.'
        )

plugin.methods.register_function(
        function=q2_moshpit.eggnog.get_references,
        inputs={},
        parameters={'mode': Str,
                    'target_taxa': List[Str],
                    'name': Str,
                    'simulate': Bool,
                    },
        outputs=[('references', ReferenceDB[Eggnog]),
                 ],
        name='get_references',
        description='Retreives necessary reference database material for'
        ' running annotations using EggnogMapper. This function will both'
        ' retrieve a complete Eggnog reference database or a specifically'
        ' selected subset of it'
        )
