# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import shutil
import numpy as np
import pandas as pd
from qiime2.util import duplicate
from .utils import _validate_mag_ids, _validate_num_partitions
from q2_types.per_sample_sequences import (
    MultiMAGSequencesDirFmt
)
from q2_types.feature_data_mag import (
    MAGSequencesDirFmt
)


def partition_feature_data_mags(
    mags: MAGSequencesDirFmt, num_partitions: int = None
) -> MAGSequencesDirFmt:
    """
    Returns a dictionary where each key is either the mag_id or an index, and
    values are the new objects with the mags.
    """
    partitioned_mags = {}
    mags_all = []

    # Get a list where every entry is a tuple representing one MAG
    for mag_id, mag_fp in mags.feature_dict().items():
        mags_all.append((mag_fp, mag_id))

    # Count number of mags and validate the num_partitions
    num_mags = len(mags_all)
    num_partitions = _validate_num_partitions(num_mags, num_partitions)
    _validate_mag_ids(num_partitions, num_mags, mags_all)

    # Split list MAGs into n arrays, where n = num_partitions
    arrays_of_mags = np.array_split(mags_all, num_partitions)

    for i, _mag in enumerate(arrays_of_mags, 1):
        result = MAGSequencesDirFmt()

        for mag_fp, mag_id in _mag:
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


def partition_sample_data_mags(
    mags: MultiMAGSequencesDirFmt, num_partitions: int = None
) -> MultiMAGSequencesDirFmt:
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
    _validate_mag_ids(num_partitions, num_mags, mags_all)

    # Split list MAGs into n arrays, where n = num_partitions
    arrays_of_mags = np.array_split(mags_all, num_partitions)

    for i, _mag in enumerate(arrays_of_mags, 1):
        result = MultiMAGSequencesDirFmt()
        samples = set([x for x, _, _ in _mag])
        print(f"Samples in this partition: {samples}")
        manifest = pd.read_csv(mags.path / "MANIFEST", index_col=None)
        print(manifest)
        manifest = manifest[manifest["sample-id"].isin(samples)]
        print(manifest)
        manifest.to_csv(result.path / "MANIFEST", index=False)

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


def collate_sample_data_mags(
    mags: MultiMAGSequencesDirFmt
) -> MultiMAGSequencesDirFmt:
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
                # Overwrite is necessary
                shutil.copy(manifest, collated_mags.path / manifest.name)

    return collated_mags


def collate_feature_data_mags(mags: MAGSequencesDirFmt) -> MAGSequencesDirFmt:
    collated_mags = MAGSequencesDirFmt()
    for mag in mags:
        for fp in mag.path.iterdir():
            duplicate(fp, collated_mags.path / fp.name)

    return collated_mags
