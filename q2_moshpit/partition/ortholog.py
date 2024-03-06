# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import shutil

from q2_types.genome_data import SeedOrthologDirFmt


def collate_orthologs(orthologs: SeedOrthologDirFmt) -> SeedOrthologDirFmt:
    result = SeedOrthologDirFmt()

    for ortholog in orthologs:
        for fp in ortholog.path.iterdir():
            shutil.move(fp, result.path / os.path.basename(fp))

    return result
