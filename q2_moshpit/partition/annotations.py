# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_types.feature_data_mag import OrthologAnnotationDirFmt
from qiime2.util import duplicate


def collate_annotations(
    ortholog_annotations: OrthologAnnotationDirFmt
) -> OrthologAnnotationDirFmt:
    # Init output
    collated_annotations = OrthologAnnotationDirFmt()

    # Copy annotations into output
    for anno in ortholog_annotations:
        for fp in anno.path.iterdir():
            duplicate(fp, collated_annotations.path / fp.name)

    return collated_annotations
