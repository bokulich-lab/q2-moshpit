# ----------------------------------------------------------------------------
# Copyright (c) COPYRIGHT_YEARS, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from q2_moshpit.plugin_setup import plugin
from q2_types.feature_data import FeatureData
from q2_types_genomics.eggnog import NOG, KEGG, OG
from q2_types_genomics.genome_data import GenomeData, Loci, Genes, Proteins

def annotate_eggnog(self, data: GenomeData[Loci | Genes | Proteins]) ->
FeatueData[NOG | OG | KEGG]:
    pass

def get_eggnog_db(self, target_database) -> None:
    pass

plugin.methods.register_function(
        function=annotate_eggnog,
        inputs={
            'data': GenomeData[Loci | Genes | Proteins],
            },
        parameters={
            },
        outputs={
            },
        name='annotate_eggnog',
        description=('Uses http://eggnog-mapper.embl.de/ to annotate'
                     'sequences')
