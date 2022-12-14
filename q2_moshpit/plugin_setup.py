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

from q2_types_genomics.per_sample_data import Contigs
from qiime2.plugin import Bool, Int
from qiime2.plugin import Plugin, Citations
# import q2_moshpit.usage_examples._examples as all_xmpls


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
            'input_sequences': 'Sequence data of the contigs we want to '
                               'search for hits using the Diamond Database',
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
        # examples={'eggnog_diamond_search':
        #           all_xmpls.eggnog_diamond_search_example},
        )

plugin.methods.register_function(
        function=q2_moshpit.annotation.eggnog_annotate_seed_orthologs,
        inputs={
            'hits_table': Ortholog[Seed],
            'eggnog_db': ReferenceDB[Eggnog],
            },
        parameters={
            'db_in_memory': Bool,
            },
        parameter_descriptions={
            'db_in_memory': 'Read entire eggnog database into memory. The '
                            'eggnog database is very large(>44GB), so this '
                            'option should only be used on clusters or other '
                            'machines with enough memory',
            },
        outputs=[('annotation_ortholog', Ortholog[Annotation])],
        name='eggnog_annotate_seed_orthologs',
        description="Uses Eggnog Mapper to apply functional annotations from "
        "the eggnog database to previously generated \"seed orthologs\".",
        # examples={'eggnog_annotate_seed_orthologs':
        #           all_xmpls.eggnog_annotate_seed_orthologs_example},
        )
