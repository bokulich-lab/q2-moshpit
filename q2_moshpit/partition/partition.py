# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from typing import Union
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


def partition_mags(
    mags: Union[MultiMAGSequencesDirFmt],
    num_partitions: int = None
) -> Union[MultiMAGSequencesDirFmt]:

    if isinstance(mags, MultiMAGSequencesDirFmt):
        return _partition_sample_data_mags(mags, num_partitions)

    elif isinstance(mags, MAGSequencesDirFmt):
        return _partition_feature_data_mags(mags, num_partitions)

    else:
        raise ValueError(
            "--i-mags is neither SampleData[MAGs] of Feature[MAGs]\n"
            "Printing associated format: \n"
            f"{type(mags)}"
        )


def collate_mags(
        mags: Union[MultiMAGSequencesDirFmt]
) -> Union[MultiMAGSequencesDirFmt]:

    if isinstance(mags[0], MultiMAGSequencesDirFmt):
        return _collate_sample_data_mags(mags)

    elif isinstance(mags[0], MAGSequencesDirFmt):
        return _collate_feature_data_mags(mags)

    else:
        raise ValueError(
            "--i-mags is neither SampleData[MAGs] of Feature[MAGs]\n"
            "Printing associated format: \n"
            f"{type(mags[0])}"
        )
