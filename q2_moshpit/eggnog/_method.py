# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import subprocess
import tempfile
from typing import Union

import pandas as pd
import qiime2.util
from typing import Union
from q2_types.per_sample_sequences import (
    ContigSequencesDirFmt, MultiMAGSequencesDirFmt
)
from q2_types.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from q2_types.reference_db import (
    EggnogRefDirFmt, DiamondDatabaseDirFmt
)
from q2_types.feature_data_mag import (
    OrthologAnnotationDirFmt, MAGSequencesDirFmt
)


def eggnog_diamond_search(
        sequences: Union[
            ContigSequencesDirFmt,
            MultiMAGSequencesDirFmt,
            MAGSequencesDirFmt,
        ],
        diamond_db: DiamondDatabaseDirFmt,
        num_cpus: int = 1, db_in_memory: bool = False
) -> (SeedOrthologDirFmt, pd.DataFrame):

    diamond_db_fp = os.path.join(str(diamond_db), 'ref_db.dmnd')
    temp = tempfile.TemporaryDirectory()

    # run analysis
    if isinstance(sequences, ContigSequencesDirFmt):
        for sample_id, contigs_fp in sequences.sample_dict().items():
            _diamond_search_runner(
                input_path=contigs_fp, diamond_db=diamond_db_fp,
                sample_label=sample_id, output_loc=temp.name,
                num_cpus=num_cpus, db_in_memory=db_in_memory
            )
    elif isinstance(sequences, MAGSequencesDirFmt):
        for mag_id, mag_fp in sequences.feature_dict().items():
            _diamond_search_runner(
                input_path=mag_fp, diamond_db=diamond_db_fp,
                sample_label=mag_id, output_loc=temp.name,
                num_cpus=num_cpus, db_in_memory=db_in_memory
            )
    elif isinstance(sequences, MultiMAGSequencesDirFmt):
        for sample_id, mags in sequences.sample_dict().items():
            for mag_id, mag_fp in mags.items():
                _diamond_search_runner(
                    input_path=mag_fp, diamond_db=diamond_db_fp,
                    sample_label=mag_id, output_loc=temp.name,
                    num_cpus=num_cpus, db_in_memory=db_in_memory
                )

    result = SeedOrthologDirFmt()
    ortholog_fps = [
        os.path.basename(x) for x
        in glob.glob(f'{temp.name}/*.seed_orthologs')
    ]
    for item in ortholog_fps:
        qiime2.util.duplicate(
            os.path.join(temp.name, item), os.path.join(result.path, item)
        )

    ft = _eggnog_feature_table(result)

    return result, ft


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


def eggnog_annotate(eggnog_hits: SeedOrthologDirFmt,
                    eggnog_db: EggnogRefDirFmt,
                    db_in_memory: bool = False,
                    num_cpus: int = 1) -> OrthologAnnotationDirFmt:

    eggnog_db_fp = eggnog_db.path

    result = OrthologAnnotationDirFmt()

    # run analysis
    for relpath, obj_path in eggnog_hits.seed_orthologs.iter_views(
            OrthologFileFmt):
        sample_label = str(relpath).rsplit(r'.', 2)[0]

        _annotate_seed_orthologs_runner(seed_ortholog=obj_path,
                                        eggnog_db=eggnog_db_fp,
                                        sample_label=sample_label,
                                        output_loc=result,
                                        db_in_memory=db_in_memory,
                                        num_cpus=num_cpus)

    return result


def _annotate_seed_orthologs_runner(seed_ortholog, eggnog_db, sample_label,
                                    output_loc, db_in_memory, num_cpus):

    # at this point instead of being able to specify the type of target
    # orthologs, we want to annotate _all_.

    cmds = ['emapper.py', '-m', 'no_search', '--annotate_hits_table',
            str(seed_ortholog), '--data_dir', str(eggnog_db),
            '-o', str(sample_label), '--output_dir', str(output_loc),
            '--cpu', str(num_cpus)]
    if db_in_memory:
        cmds.append('--dbmem')

    subprocess.run(cmds, check=True)
