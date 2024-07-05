# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess

from q2_types.feature_data import BLAST6
from q2_types.feature_data_mag import OrthologAnnotationDirFmt
from q2_types.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from q2_types.reference_db import EggnogRefDirFmt
from q2_types.sample_data import SampleData


def _annotate_seed_orthologs_runner(
        seed_ortholog, eggnog_db, sample_label, output_loc,
        db_in_memory, num_cpus
):
    # at this point instead of being able to specify the type of target
    # orthologs, we want to annotate _all_.

    cmds = ['emapper.py', '-m', 'no_search', '--annotate_hits_table',
            str(seed_ortholog), '--data_dir', str(eggnog_db),
            '-o', str(sample_label), '--output_dir', str(output_loc),
            '--cpu', str(num_cpus)]
    if db_in_memory:
        cmds.append('--dbmem')

    subprocess.run(cmds, check=True, cwd=output_loc)


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
    collate_annotations = ctx.get_action("moshpit", "collate_annotations")

    if eggnog_hits.type <= SampleData[BLAST6]:
        partition_method = ctx.get_action("moshpit", "partition_orthologs")
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
