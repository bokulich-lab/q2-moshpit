# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from q2_types.sample_data import SampleData
from  q2_types_genomics.per_sample_data import ContigSequencesDirFmt, Contigs
from q2_types_genomics.reference_db import ReferenceDB, Diamond, Eggnog
from q2_types_genomics.eggnog import EggnogRefDirFmt
from q2_types_genomics.ortholog import OrthologDirFmt, Ortholog, Seed
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt, Contigs
from q2_types.feature_data import DNAFASTAFormat, FeatureData
from q2_types_genomics.diamond import DiamondDatabaseDirFmt
from q2_moshpit.plugin_setup import plugin

import os
import subprocess
import shutil

# diamond search
def eggnog_search_diamond(input_sequences: ContigSequencesDirFmt,
                   diamond_db:DiamondDatabaseDirFmt,
                   eggnog_db:EggnogRefDirFmt,
                   )-> OrthologDirFmt:

    result = OrthologDirFmt()
    diamond_fp = diamond_db.path / 'ref_db.dmnd'

    for relpath, obj_path in input_sequences.sequences.iter_views(DNAFASTAFormat):
        sample_label = str(relpath).rsplit(r'\.', 2)[0] #+ '_seed_ortholog'

        _diamond_search_runner(input_path=obj_path,
                               eggnog_db=str(eggnog_db),
                               diamond_db=diamond_fp,
                               sample_label=sample_label,
                               output=str(result))

    return result

def _diamond_search_runner(input_path, eggnog_db, diamond_db, sample_label, output):
    cmds = ['emapper.py',
            '--data_dir', eggnog_db,
            '--dmnd_db', diamond_db,
            '-i', str(input_path),
            '--itype', 'metagenome',
            '-o', sample_label,
            '--output_dir', output,
            '--no_annot',
            ]
    subprocess.run(cmds, check=True)
    return None

