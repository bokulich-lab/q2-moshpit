# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from typing import Union
from .utils import (
    _collate_sample_data_mags, _collate_feature_data_mags,
    _partition_feature_data_mags, _partition_sample_data_mags
)
from q2_types_genomics.per_sample_data._format import (
    MultiMAGSequencesDirFmt
)
from q2_types_genomics.feature_data._format import (
    MAGSequencesDirFmt
)


def partition_mags(
    mags: Union(MultiMAGSequencesDirFmt, MAGSequencesDirFmt),
    num_partitions: int = None
) -> Union(MultiMAGSequencesDirFmt, MAGSequencesDirFmt):

    if isinstance(mags, MultiMAGSequencesDirFmt):
        _partition_sample_data_mags(mags, num_partitions)

    elif isinstance(mags, MAGSequencesDirFmt):
        _partition_feature_data_mags(mags, num_partitions)

    else:
        raise ValueError(
            "--i-mags is neither SampleData[MAGs] of Feature[MAGs]\n"
            "Printing associated format: \n"
            f"{type(mags)}"
        )


def collate_mags(
        mags: Union(MultiMAGSequencesDirFmt, MAGSequencesDirFmt)
) -> Union(MultiMAGSequencesDirFmt, MAGSequencesDirFmt):

    if isinstance(mags, MultiMAGSequencesDirFmt):
        return _collate_sample_data_mags(mags)

    elif isinstance(mags, MAGSequencesDirFmt):
        return _collate_feature_data_mags(mags)

    else:
        raise ValueError(
            "--i-mags is neither SampleData[MAGs] of Feature[MAGs]\n"
            "Printing associated format: \n"
            f"{type(mags)}"
        )
