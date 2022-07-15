# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
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


import subprocess
def run_commands(cmds, verbose=True):
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)


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
        )
