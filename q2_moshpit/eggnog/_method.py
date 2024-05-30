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

from q2_types.feature_data import FeatureData, BLAST6
from q2_types.feature_data_mag import (
    OrthologAnnotationDirFmt, MAGSequencesDirFmt, MAG
)
from q2_types.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from q2_types.per_sample_sequences import (
    ContigSequencesDirFmt, MultiMAGSequencesDirFmt, Contigs, MAGs
)
from q2_types.reference_db import (
    EggnogRefDirFmt, DiamondDatabaseDirFmt
)
from q2_types.sample_data import SampleData


def eggnog_diamond_search(
        ctx,
        sequences,
        diamond_db,
        num_cpus=1,
        db_in_memory=False,
        num_partitions=None
):
    _eggnog_diamond_search = ctx.get_action(
        "moshpit", "_eggnog_diamond_search")

    if sequences.type <= FeatureData[MAG]:
        partition_method = ctx.get_action("types", "partition_feature_data_mags")
    elif sequences.type <= SampleData[Contigs]:
        partition_method = ctx.get_action("assembly", "partition_contigs")
    elif sequences.type <= SampleData[MAGs]:
        partition_method = ctx.get_action("types", "partition_sample_data_mags")
    else:
        raise NotImplementedError()

    collate_hits = ctx.get_action("types", "collate_ortholog_annotations")
    _eggnog_feature_table = ctx.get_action("moshpit", "_eggnog_feature_table")

    (partitioned_sequences,) = partition_method(sequences, num_partitions)

    hits = []
    for seq in partitioned_sequences.values():
        (hit, _) = _eggnog_diamond_search(
            seq, diamond_db, num_cpus, db_in_memory)
        hits.append(hit)

    (collated_hits,) = collate_hits(hits)
    (collated_tables,) = _eggnog_feature_table(collated_hits)

    return collated_hits, collated_tables


def _eggnog_diamond_search(
        sequences: Union[
            ContigSequencesDirFmt,
            MultiMAGSequencesDirFmt,
            MAGSequencesDirFmt
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


def _eggnog_annotate(
        eggnog_hits: SeedOrthologDirFmt,
        eggnog_db: EggnogRefDirFmt,
        db_in_memory: bool = False,
        num_cpus: int = 1
) -> OrthologAnnotationDirFmt:

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


def eggnog_annotate(
        ctx,
        eggnog_hits,
        eggnog_db,
        db_in_memory=False,
        num_cpus=1,
        num_partitions=None
):
    _eggnog_annotate = ctx.get_action("moshpit", "_eggnog_annotate")
    collate_annotations = ctx.get_action("types", "collate_annotations")

    if eggnog_hits.type <= SampleData[BLAST6]:
        partition_method = ctx.get_action("types", "partition_orthologs")
    else:
        raise NotImplementedError()

    (partitioned_orthologs, ) = partition_method(eggnog_hits, num_partitions)

    annotations = []
    for orthologs in partitioned_orthologs.values():
        (annotation, ) = _eggnog_annotate(
            orthologs, eggnog_db, db_in_memory, num_cpus
        )
        annotations.append(annotation)

    (collated_annotations, ) = collate_annotations(annotations)
    return collated_annotations


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
