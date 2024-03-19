# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from typing import List, Dict, Tuple
import warnings


def _validate_num_partitions(
        num_samples: int, num_partitions: int, sample_type: str = "sample"
) -> int:

    if num_partitions is None:
        return num_samples
    elif num_partitions > num_samples:
        warnings.warn(
            "You have requested a number of partitions "
            f"'{num_partitions}' that is greater than your number "
            f"of {sample_type}s '{num_samples}.' Your data will be "
            f"partitioned by {sample_type} into '{num_samples}' "
            "partitions."
        )
        return num_samples
    else:
        return num_partitions


def _validate_mag_ids(
    num_partitions: int, num_mags: int, mags_all: List[tuple]
):
    # If num_partitions == num_mags and MAG ids are not unique
    # the output will be missing these duplicated-id MAGs.
    # While this is technically impossible since
    # MAGs should have unique IDs by construction, it could still happen that a
    # used imports MAGs with non-unique IDs. In such case this test would be
    # useful.

    if num_partitions == num_mags:
        mag_ids = [mag_id[1] for mag_id in mags_all]
        duplicates = [
            mag_id for mag_id in mag_ids if mag_ids.count(mag_id) > 1
        ]
        if len(duplicates) > 0:
            raise ValueError(
                "MAG IDs are not unique. "
                "They must be unique in order to output all partitions "
                "correctly. Printing duplicate MAG IDs: "
                f"{set(duplicates)}"
            )
