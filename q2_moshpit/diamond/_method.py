# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from q2_types_genomics.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from q2_types.feature_data import DNAFASTAFormat
from q2_types_genomics.reference_db import DiamondDatabaseDirFmt
import qiime2.util

import os
import subprocess
import shutil
import tempfile
import re
import pandas as pd


def eggnog_diamond_search(input_sequences: ContigSequencesDirFmt,
                          diamond_db: DiamondDatabaseDirFmt,
                          num_cpus: int = 1, db_in_memory: bool = False
                          ) -> (SeedOrthologDirFmt, pd.DataFrame):

    diamond_db_fp = os.path.join(str(diamond_db), 'ref_db.dmnd')
    temp = tempfile.TemporaryDirectory()

    # run analysis
    for relpath, obj_path in input_sequences.sequences.iter_views(
            DNAFASTAFormat):
        sample_label = str(relpath).rsplit(r'_', 1)[0]

        _diamond_search_runner(input_path=obj_path,
                               diamond_db=diamond_db_fp,
                               sample_label=sample_label,
                               output_loc=temp.name,
                               num_cpus=num_cpus,
                               db_in_memory=db_in_memory)

    result = SeedOrthologDirFmt()

    for item in os.listdir(temp.name):
        if re.match(r".*\.seed_orthologs", item):
            qiime2.util.duplicate(os.path.join(temp.name, item),
                                  os.path.join(result.path, item))

    ft = extract_ft_from_seed_orthologs(result)

    return (result, ft)


def extract_ft_from_seed_orthologs(seed_orthologs: SeedOrthologDirFmt
                                   ) -> pd.DataFrame:

    per_sample_counts = []

    for sample_path, obj in seed_orthologs.seed_orthologs.iter_views(
            OrthologFileFmt):
        #TODO: put filename to sample name logic on OrthologFileFmt object
        sample_name = str(sample_path).replace('.emapper.seed_orthologs', '')
        sample_df = obj.view(pd.DataFrame)
        sample_feature_counts = sample_df.value_counts('sseqid')
        sample_feature_counts.name = str(sample_name)
        per_sample_counts.append(sample_feature_counts)
    df = pd.DataFrame(per_sample_counts)
    df.fillna(0, inplace=True)
    df.columns = df.columns.astype('str')

    return df


def _diamond_search_runner(input_path, diamond_db, sample_label, output_loc,
                           num_cpus, db_in_memory):

    cmds = ['emapper.py', '-i', str(input_path), '-o', sample_label,
            '-m', 'diamond', '--no_annot', '--dmnd_db', str(diamond_db),
            '--itype', 'metagenome', '--output_dir', output_loc, '--cpu',
            str(num_cpus)]
    if db_in_memory:
        cmds.append('--dbmem')

    subprocess.run(cmds, check=True)
