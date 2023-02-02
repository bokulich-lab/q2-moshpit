# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


# plugin imports
from q2_types_genomics.genome_data import OrthologFileFmt, SeedOrthologDirFmt
from q2_types_genomics.feature_data import OrthologAnnotationDirFmt
from q2_types_genomics.reference_db import EggnogRefDirFmt
from q2_types_genomics.reference_db import ReferenceDB, EggnogRefDirFmt

# library imports
import subprocess
import shutil
import tempfile
import os


def eggnog_annotate_seed_orthologs(hits_table: SeedOrthologDirFmt,
                                   eggnog_db: EggnogRefDirFmt,
                                   db_in_memory: bool = False,
                                   ) -> OrthologAnnotationDirFmt:

    eggnog_db_fp = eggnog_db.path
    temp = tempfile.TemporaryDirectory()

    # run analysis
    for relpath, obj_path in hits_table.seed_orthologs.iter_views(
            OrthologFileFmt):
        sample_label = str(relpath).rsplit(r'\.', 2)[0] + '.annotation_ortholog'

        _annotate_seed_orthologs_runner(seed_ortholog=obj_path,
                                        eggnog_db=eggnog_db_fp,
                                        sample_label=sample_label,
                                        output_loc=temp.name,
                                        db_in_memory=db_in_memory)

    # INSTANTIATE RESULT OBJECT
    result = OrthologAnnotationDirFmt()

    for item in os.listdir(temp.name):
        shutil.copy(os.path.join(temp.name, item), result.path)

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
