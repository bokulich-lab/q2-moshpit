# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


#plugin imports
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from q2_types_genomics.ortholog import AnnotationOrthologDirFmt, SeedOrthologDirFmt
from q2_types_genomics.eggnog import EggnogRefDirFmt

#library imports
import os
import subprocess
import shutil
import tempfile
import re


def eggnog_annotate_seed_orthologs(hits_table: SeedOrthologDirFmt,
                                   eggnog_db: EggnogRefDirFmt,
                                   ) -> AnnotationOrthologDirFmt:

    diamond_db_fp = diamond_db.path / 'ref_db.dmnd'
    temp = tempfile.TemporaryDirectory()

    # run analysis
    for relpath, obj_path in input_sequences.sequences.iter_views(
            DNAFASTAFormat):
        sample_label = str(relpath).rsplit(r'\.', 2)[0] + '_seed_ortholog'

        # RUN HELPER
        # _diamond_search_runner(input_path=obj_path,
        #                        diamond_db=diamond_db_fp,
        #                        sample_label=sample_label,
        #                        output_loc=temp.name)

    # INSTANTIATE RESULT OBJECT
    result = AnnotationOrthologDirFmt()

    # COPY TO RESULT ARTIFACT
    # for item in os.listdir(temp.name):
    #     if re.match(r".*\..*\.seed_orthologs", item):
    #         shutil.copy(os.path.join(temp.name, item), result.path)

    return result


def _annotate_seed_orthologs_runner(input_path, eggnog_db, sample_label, output_loc):

    # at this point instead of being able to specify the type of target
    # orthologs, we want to annotate _all_.

    cmds = ['emapper.py', '-m', 'no_search', '--annotate_hits_table', str(hits_table), '--data_dir', str(eggnog_db),
            '-o', str(sample_label), '--output_dir', str(output_loc)]


    subprocess.run(cmds, check=True)
