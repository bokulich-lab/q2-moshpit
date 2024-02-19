# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from .partition import (
    collate_feature_data_mags, collate_sample_data_mags,
    partition_feature_data_mags, partition_sample_data_mags
)

__all__ = [
    "collate_feature_data_mags", "collate_sample_data_mags",
    "partition_feature_data_mags", "partition_sample_data_mags"
]
