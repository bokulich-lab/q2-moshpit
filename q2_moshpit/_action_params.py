# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.core.type import Bool, Int, Range, Float, Str

checkm_params = {
    'db_path': Str,
    'reduced_tree': Bool,
    'unique': Int % Range(1, None),
    'multi': Int % Range(1, None),
    'force_domain': Bool,
    'no_refinement': Bool,
    'individual_markers': Bool,
    'skip_adj_correction': Bool,
    'skip_pseudogene_correction': Bool,
    'aai_strain': Float % Range(0, 1),
    'ignore_thresholds': Bool,
    'e_value': Float % Range(0, 1),
    'length': Float % Range(0, 1),
    'threads': Int % Range(1, None),
    'pplacer_threads': Int % Range(1, None)
}

checkm_param_descriptions = {
    'db_path': 'Path to the database required by CheckM. For more '
               'details see: https://github.com/Ecogenomics/CheckM/'
               'wiki/Installation#required-reference-data.',
    'reduced_tree': 'Use reduced tree (requires <16GB of memory) for '
                    'determining lineage of each bin.',
    'unique': 'Minimum number of unique phylogenetic markers required to use '
              'lineage-specific marker set. Default: 10.',
    'multi': 'Maximum number of multi-copy phylogenetic markers before '
             'defaulting to domain-level marker set. Default: 10.',
    'force_domain': 'Use domain-level sets for all bins.',
    'no_refinement': 'Do not perform lineage-specific marker set refinement.',
    'individual_markers': 'Treat marker as independent '
                          '(i.e., ignore co-located set structure).',
    'skip_adj_correction': 'Do not exclude adjacent marker genes when '
                           'estimating contamination.',
    'skip_pseudogene_correction': 'Skip identification and filtering '
                                  'of pseudogenes.',
    'aai_strain': 'AAI threshold used to identify strain heterogeneity. '
                  'Default: 0.9.',
    'ignore_thresholds': 'Ignore model-specific score thresholds.',
    'e_value': 'E-value cut off. Default: 1e-10.',
    'length': 'Percent overlap between target and query. Default: 0.7.',
    'threads': 'Number of threads. Default: 1.',
    'pplacer_threads': 'Number of threads used by pplacer (memory usage '
                       'increases linearly with additional threads). '
                       'Default: 1.'
}
