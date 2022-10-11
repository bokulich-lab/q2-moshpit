# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_types_genomics.eggnog import Ortholog, Seed, ArbitraryHeaderTSVDirFmt, ReferenceDB, Eggnog, Diamond
from q2_types.feature_data import DNAFASTAFormat, FeatureData

from q2_types.sample_data import SampleData
from q2_types_genomics.per_sample_data import Contigs

from ..plugin_setup import plugin

# diamond search
# inputs DNAFASTAFormat data, eggnog reference database, diamond reference database
# outputs: seed orthologs(fmt: ArbitraryHeaderTSVDirFmt, type: Ortholog[Seed])

def diamond_search(input_sequences: DNAFASTAFormat, diamond_db, eggnog_db, seed_orthologs=None):
    pass
