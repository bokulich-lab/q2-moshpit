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

import pandas as pd
import qiime2.util

from q2_types.feature_data import FeatureData
from q2_types.feature_data_mag import MAG, MAGSequencesDirFmt
from q2_types.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from q2_types.per_sample_sequences import (
    Contigs, MAGs, ContigSequencesDirFmt, MultiMAGSequencesDirFmt
)
from q2_types.sample_data import SampleData


def _symlink_files_to_target_dir(
    pressed_hmm_db, idmap, seed_alignments, target_dir
):
    for source_dir in [str(pressed_hmm_db), str(idmap), str(seed_alignments)]:
        for filename in os.listdir(source_dir):
            os.symlink(
                os.path.join(source_dir, filename),
                os.path.join(target_dir, filename)
            )


def _run_eggnog_search_pipeline(
    ctx, sequences, db, num_cpus, db_in_memory, num_partitions, search_action
):
    if sequences.type <= FeatureData[MAG]:
        plugin, action_name = "moshpit", "partition_feature_data_mags"
    elif sequences.type <= SampleData[Contigs]:
        plugin, action_name = "assembly", "partition_contigs"
    elif sequences.type <= SampleData[MAGs]:
        plugin, action_name = "moshpit", "partition_sample_data_mags"
    else:
        raise NotImplementedError()

    partition_method = ctx.get_action(plugin, action_name)
    _eggnog_search = ctx.get_action("moshpit", search_action)
    collate_hits = ctx.get_action("moshpit", "collate_orthologs")
    _eggnog_feature_table = ctx.get_action("moshpit", "_eggnog_feature_table")
    (partitioned_sequences,) = partition_method(sequences, num_partitions)

    hits = []
    for seq in partitioned_sequences.values():
        (hit, _) = _eggnog_search(seq, *db, num_cpus, db_in_memory)
        hits.append(hit)

    (collated_hits,) = collate_hits(hits)
    (collated_tables,) = _eggnog_feature_table(collated_hits)
    return collated_hits, collated_tables


def _search_runner(
    input_path, sample_label, output_loc, num_cpus, db_in_memory,
    runner_args
):
    cmd = [
        'emapper.py', '-i', str(input_path), '-o', sample_label,
        '-m', *runner_args,
        '--itype', 'metagenome', '--output_dir', output_loc,
        '--cpu', str(num_cpus), '--no_annot'
    ]
    if db_in_memory:
        cmd.append('--dbmem')

    try:
        subprocess.run(cmd, check=True, cwd=output_loc)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "Error running eggNOG-mapper. "
            f"The exception was: {e}"
        )


def _eggnog_search(
    sequences, search_runner, output_loc
) -> (SeedOrthologDirFmt, pd.DataFrame):  # type: ignore
    # run analysis
    if isinstance(sequences, ContigSequencesDirFmt):
        for sample_id, contigs_fp in sequences.sample_dict().items():
            search_runner(input_path=contigs_fp, sample_label=sample_id)
    elif isinstance(sequences, MAGSequencesDirFmt):
        for mag_id, mag_fp in sequences.feature_dict().items():
            search_runner(input_path=mag_fp, sample_label=mag_id)
    elif isinstance(sequences, MultiMAGSequencesDirFmt):
        for sample_id, mags in sequences.sample_dict().items():
            for mag_id, mag_fp in mags.items():
                search_runner(input_path=mag_fp, sample_label=mag_id)

    result = SeedOrthologDirFmt()
    ortholog_fps = [
        os.path.basename(x) for x
        in glob.glob(f'{output_loc}/*.seed_orthologs')
    ]
    for item in ortholog_fps:
        qiime2.util.duplicate(
            os.path.join(output_loc, item),
            os.path.join(result.path, item)
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
