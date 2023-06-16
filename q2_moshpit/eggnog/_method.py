# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import subprocess
import os
import tempfile
import re
import pandas as pd

from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from q2_types_genomics.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from q2_types_genomics.feature_data import OrthologAnnotationDirFmt
from q2_types_genomics.reference_db import EggnogRefDirFmt
from q2_types.feature_data import DNAFASTAFormat
from q2_types_genomics.reference_db import DiamondDatabaseDirFmt
import qiime2.util


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

    ft = _eggnog_feature_table(result)

    return (result, ft)


def _eggnog_feature_table(seed_orthologs: SeedOrthologDirFmt) -> pd.DataFrame:

    per_sample_counts = []

    for sample_path, obj in seed_orthologs.seed_orthologs.iter_views(
            OrthologFileFmt):
        # TODO: put filename to sample name logic on OrthologFileFmt object
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


def eggnog_annotate(hits_table: SeedOrthologDirFmt,
                    eggnog_db: EggnogRefDirFmt,
                    db_in_memory: bool = False) -> OrthologAnnotationDirFmt:

    eggnog_db_fp = eggnog_db.path

    result = OrthologAnnotationDirFmt()

    # run analysis
    for relpath, obj_path in hits_table.seed_orthologs.iter_views(
            OrthologFileFmt):
        sample_label = str(relpath).rsplit(r'.', 2)[0]

        _annotate_seed_orthologs_runner(seed_ortholog=obj_path,
                                        eggnog_db=eggnog_db_fp,
                                        sample_label=sample_label,
                                        output_loc=result,
                                        db_in_memory=db_in_memory)

    return result


def _annotate_seed_orthologs_runner(seed_ortholog, eggnog_db, sample_label,
                                    output_loc, db_in_memory):

    # at this point instead of being able to specify the type of target
    # orthologs, we want to annotate _all_.

    cmds = ['emapper.py', '-m', 'no_search', '--annotate_hits_table',
            str(seed_ortholog), '--data_dir', str(eggnog_db),
            '-o', str(sample_label), '--output_dir', str(output_loc)]
    if db_in_memory:
        cmds.append('--dbmem')

    subprocess.run(cmds, check=True)
