# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import warnings

import numpy as np
from qiime2.util import duplicate

from q2_types_genomics.per_sample_data._format import (
    MultiFASTADirectoryFormat
)


def partition_mags(
    mags: MultiFASTADirectoryFormat, num_partitions: int = None
) -> MultiFASTADirectoryFormat:
    partitioned_mags = {}
    mags_all = []
    for sample_id, mag in mags.sample_dict().items():
        for mag_id, mag_fp in mag.items():
            mags_all.append((sample_id, mag_id, mag_fp))

    num_mags = len(mags_all)
    if num_partitions is None:
        num_partitions = num_mags
    elif num_partitions > num_mags:
        warnings.warn(
            "You have requested a number of partitions"
            f" '{num_partitions}' that is greater than your number"
            f" of MAGs '{num_mags}.' Your data will be"
            f" partitioned by MAG into '{num_mags}'"
            " partitions."
        )
        num_partitions = num_mags

    mags = np.array_split(mags_all, num_partitions)
    for i, _mag in enumerate(mags, 1):
        result = MultiFASTADirectoryFormat()

        for sample_id, mag_id, mag_fp in _mag:
            os.makedirs(result.path / sample_id, exist_ok=True)
            duplicate(
                mag_fp,
                result.path / sample_id / os.path.basename(mag_fp)
            )

        # If num_partitions == num_mags we will only have gone through one
        # MAG in the above loop and will use its id as a key. Otherwise we
        # may have gone through multiple MAGs in the above loop and will be
        # using indices for keys
        if num_partitions == num_mags:
            partitioned_mags[mag_id] = result
        else:
            partitioned_mags[i] = result

    return partitioned_mags


def collate_mags(
        mags: MultiFASTADirectoryFormat
) -> MultiFASTADirectoryFormat:
    collated_mags = MultiFASTADirectoryFormat()

    for mag in mags:
        for fp in mag.path.iterdir():
            duplicate(fp, collated_mags.path / fp.name)

    return collated_mags
