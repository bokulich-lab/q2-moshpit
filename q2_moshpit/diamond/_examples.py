# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os

def cook_contigs_factory():
    import qiime2
    return qiime2.Artifact.load(contig_fp)

def cook_diamond_db_factory():
    import qiime2
    return qiime2.Artifact.load(diamond_db_fp)

def eggnog_search_diamond(use):
    contigs = use.init_artifact('contigs', cook_contigs_factory)
    diamond_db = use.init_artifact('diamond_db', cook_diamond_db_factory)

    seed_ortholog_return = use.action(
            use.UsageAction(plugin_id='moshpit', action_id='eggnog-search-diamond'),
            use.UsageInputs(input_sequences=contigs, diamond_db=diamond_db),
            use.UsageOutputNames(seed_ortholog='seed_ortholog')
            )
    # remove after failure, this is just a guard to ensure running
    seed_ortholog_return.assert_output_type('Visualization')
    seed_ortholog_return.assert_output_type('Artifact')
