# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from q2_types_genomics.ortholog import SeedOrthologDirFmt
from q2_types.feature_data import DNAFASTAFormat
from q2_types_genomics.diamond import DiamondDatabaseDirFmt

import os
import subprocess
import shutil
import tempfile
import re


def eggnog_diamond_search(input_sequences: ContigSequencesDirFmt,
                          diamond_db: DiamondDatabaseDirFmt,
                          num_cpus: int = 1,
                          ) -> SeedOrthologDirFmt:

    diamond_db_fp = os.path.join(str(diamond_db), 'ref_db.dmnd')
    temp = tempfile.TemporaryDirectory()

    # run analysis
    for relpath, obj_path in input_sequences.sequences.iter_views(
            DNAFASTAFormat):
        sample_label = str(relpath).rsplit(r'\.', 2)[0] + '_seed_ortholog'

        _diamond_search_runner(input_path=obj_path,
                               diamond_db=diamond_db_fp,
                               sample_label=sample_label,
                               output_loc=temp.name,
                               num_cpus=num_cpus)

    result = SeedOrthologDirFmt()

    for item in os.listdir(temp.name):
        if re.match(r".*\..*\.seed_orthologs", item):
            shutil.copy(os.path.join(temp.name, item), result.path)

    return result


def _diamond_search_runner(input_path, diamond_db, sample_label, output_loc, num_cpus):

    cmds = ['emapper.py', '-i', str(input_path), '-o', sample_label,
            '-m', 'diamond', '--no_annot', '--dmnd_db', str(diamond_db),
            '--itype', 'metagenome', '--output_dir', output_loc, '--cpu', str(num_cpus)]

    subprocess.run(cmds, check=True)
