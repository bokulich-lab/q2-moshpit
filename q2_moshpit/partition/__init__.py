# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from .annotations import collate_annotations
from .mags import (
    collate_feature_data_mags, collate_sample_data_mags,
    partition_feature_data_mags, partition_sample_data_mags,

)
from .ortholog import (
    collate_orthologs, partition_orthologs
)

__all__ = [
    "collate_feature_data_mags", "collate_sample_data_mags",
    "partition_feature_data_mags", "partition_sample_data_mags",
    "collate_orthologs", "partition_orthologs", "collate_annotations"
]
