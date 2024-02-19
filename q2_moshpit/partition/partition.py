# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_moshpit.partition.utils import (
    _collate_sample_data_mags, _collate_feature_data_mags,
    _partition_feature_data_mags, _partition_sample_data_mags
)
from q2_types_genomics.per_sample_data import (
    MultiMAGSequencesDirFmt
)
from q2_types_genomics.feature_data import (
    MAGSequencesDirFmt
)


def partition_feature_data_mags(
    mags: MAGSequencesDirFmt, num_partitions: int = None
) -> MAGSequencesDirFmt:
    return _partition_feature_data_mags(mags, num_partitions)


def partition_sample_data_mags(
    mags: MultiMAGSequencesDirFmt, num_partitions: int = None
) -> MultiMAGSequencesDirFmt:
    return _partition_sample_data_mags(mags, num_partitions)


def collate_sample_data_mags(
    mags: MultiMAGSequencesDirFmt
) -> MultiMAGSequencesDirFmt:
    return _collate_sample_data_mags(mags)


def collate_feature_data_mags(mags: MAGSequencesDirFmt) -> MAGSequencesDirFmt:
    return _collate_feature_data_mags(mags)
