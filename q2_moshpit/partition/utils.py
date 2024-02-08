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
    MultiMAGSequencesDirFmt
)
from q2_types_genomics.feature_data._format import (
    MAGSequencesDirFmt
)


def _validate_num_partitions(num_mags, num_partitions) -> int:
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
    return num_partitions


def _partition_feature_data_mags(mags, num_partitions) -> dict:
    """
    Returns a dictionary where each key is either the mag_id or an index, and
    values are the new objects with the mags.
    """
    partitioned_mags = {}
    mags_all = []

    # Get a list where every entry is a tuple representing one MAG
    for mag_id, mag_fp in mags.feature_dict().items():
        mags_all.append((mag_id, mag_fp))

    # Count number of mags and validate the num_partitions
    num_mags = len(mags_all)
    num_partitions = _validate_num_partitions(num_mags, num_partitions)

    # Split list MAGs into n arrays, where n = num_partitions
    arrays_of_mags = np.array_split(mags_all, num_partitions)

    for i, _mag in enumerate(arrays_of_mags, 1):
        result = MAGSequencesDirFmt()

        for mag_id, mag_fp in _mag:
            duplicate(mag_fp, result.path / os.path.basename(mag_fp))

        # If num_partitions == num_mags we will only have gone through one
        # MAG in the above loop and will use its id as a key. Otherwise we
        # may have gone through multiple MAGs in the above loop and will be
        # using indices for keys
        if num_partitions == num_mags:
            partitioned_mags[mag_id] = result
        else:
            partitioned_mags[i] = result

    return partitioned_mags


def _partition_sample_data_mags(mags, num_partitions) -> dict:
    """
    Returns a dictionary where each key is either the mag_id or an index, and
    values are the new objects with the mags.
    """
    partitioned_mags = {}
    mags_all = []

    # Get a list where every entry is a tuple representing one MAG
    for sample_id, mag in mags.sample_dict().items():
        for mag_id, mag_fp in mag.items():
            mags_all.append((sample_id, mag_id, mag_fp))

    # Count number of mags and validate the num_partitions
    num_mags = len(mags_all)
    num_partitions = _validate_num_partitions(num_mags, num_partitions)

    # Split list MAGs into n arrays, where n = num_partitions
    arrays_of_mags = np.array_split(mags_all, num_partitions)

    for i, _mag in enumerate(arrays_of_mags, 1):
        result = MultiMAGSequencesDirFmt()
        duplicate(mags.path / "MANIFEST", result.path / "MANIFEST")

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


def _collate_sample_data_mags(mags):
    collated_mags = MultiMAGSequencesDirFmt()

    # For every partition
    for mag in mags:

        # For every sample in the partition
        for file_or_dir in mag.path.iterdir():

            if file_or_dir.is_dir():
                sample = file_or_dir
                os.makedirs(collated_mags.path / sample.name, exist_ok=True)

                # For every mag in the sample
                for mag in sample.iterdir():
                    duplicate(mag, collated_mags.path / sample.name / mag.name)

            # If its a file, it should be the manifest
            # Since its present many times it will be overwritten, but that ok
            else:
                manifest = file_or_dir
                duplicate(manifest, collated_mags.path / manifest.name)

    return collated_mags


def _collate_feature_data_mags(mags):
    collated_mags = MAGSequencesDirFmt()
    for mag in mags:
        for fp in mag.path.iterdir():
            duplicate(fp, collated_mags.path / fp.name)

    return collated_mags
