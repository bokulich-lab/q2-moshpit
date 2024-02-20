# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import filecmp
from unittest.mock import patch
from q2_moshpit.partition.partition import (
    partition_sample_data_mags, partition_feature_data_mags,
    collate_feature_data_mags, collate_sample_data_mags
)
from qiime2.plugin.testing import TestPluginBase
from q2_types_genomics.per_sample_data._format import (
    MultiMAGSequencesDirFmt,
)
from q2_types_genomics.feature_data._format import (
    MAGSequencesDirFmt
)


class TestHelpers(TestPluginBase):
    package = 'q2_moshpit.partition.tests'

    @patch('q2_moshpit.partition.utils._validate_mag_ids')
    @patch('q2_moshpit.partition.utils._validate_num_partitions')
    def test_partition_sample_data_mags(
        self,
        mock_validate_num_partitions,
        mock_validate_mag_ids
    ):
        # Partition mags
        p = self.get_data_path("collated/sample_data")
        mags = MultiMAGSequencesDirFmt(path=p, mode="r")
        mock_validate_num_partitions.return_value = 2
        partitioned_mags = partition_sample_data_mags(mags, 2)

        # Expected mag ids for every sample
        mag_ids_sample_1 = [
            "24dee6fe-9b84-45bb-8145-de7b092533a1.fasta",
            "fb0bc871-04f6-486b-a10e-8e0cb66f8de3.fasta"
        ]
        mag_ids_sample_2 = [
            "d65a71fa-4279-4588-b937-0747ed5d604d.fasta",
        ]

        # Compare dirs
        for i, mag_ids in [(1, mag_ids_sample_1), (2, mag_ids_sample_2)]:
            dircmp = filecmp.dircmp(
                partitioned_mags[i].path,
                mags.path
            )
            self.assertListEqual(
                ["MANIFEST", f"sample{i}"], dircmp.common
            )
            dircmp = filecmp.dircmp(
                f"{partitioned_mags[i].path}/sample{i}",
                f"{mags.path}/sample{i}"
            )
            self.assertListEqual(
                [
                    *mag_ids,
                ],
                dircmp.common
            )

    @patch('q2_moshpit.partition.utils._validate_mag_ids')
    @patch('q2_moshpit.partition.utils._validate_num_partitions')
    def test_partition_feature_data_mags(
        self,
        mock_validate_num_partitions,
        mock_validate_mag_ids
    ):
        # Partition Feature Data
        p = self.get_data_path("collated/feature_data")
        mags = MAGSequencesDirFmt(path=p, mode="r")
        mock_validate_num_partitions.return_value = 2
        partitioned_mags = partition_feature_data_mags(mags)

        # Expected mag ids
        mag_ids = [
            "24dee6fe-9b84-45bb-8145-de7b092533a1",
            "fb0bc871-04f6-486b-a10e-8e0cb66f8de3"
        ]

        # compare partitions
        for i in [0, 1]:
            dircmp = filecmp.dircmp(
                partitioned_mags[mag_ids[i]].path, mags.path
            )
            self.assertListEqual([f"{mag_ids[i]}.fasta"], dircmp.common)

    def test_collate_sample_data_mags(self):
        p1 = self.get_data_path("partitioned/sample_data/mag1")
        p2 = self.get_data_path("partitioned/sample_data/mag2")
        mags = [
            MultiMAGSequencesDirFmt(p1, mode="r"),
            MultiMAGSequencesDirFmt(p2, mode="r")
        ]

        collated_mags = collate_sample_data_mags(mags)
        expected = self.get_data_path("collated/sample_data")

        # compare first dir
        dircmp = filecmp.dircmp(collated_mags.path, expected)
        self.assertListEqual(["MANIFEST", "sample1", "sample2"], dircmp.common)

        # Compare second dir
        dircmp = filecmp.dircmp(
            f"{collated_mags.path}/sample1",
            f"{expected}/sample1"
        )
        self.assertListEqual(
            [
                "24dee6fe-9b84-45bb-8145-de7b092533a1.fasta",
                "fb0bc871-04f6-486b-a10e-8e0cb66f8de3.fasta"
            ],
            dircmp.common
        )

        # compare third dir
        dircmp = filecmp.dircmp(
            f"{collated_mags.path}/sample2",
            f"{expected}/sample2"
        )
        self.assertListEqual(
            ["d65a71fa-4279-4588-b937-0747ed5d604d.fasta"],
            dircmp.common
        )

    def test_collate_feature_data_mags(self):
        # collate test data
        p1 = self.get_data_path("partitioned/feature_data/mag1")
        p2 = self.get_data_path("partitioned/feature_data/mag2")
        mags = [
            MAGSequencesDirFmt(p1, mode="r"),
            MAGSequencesDirFmt(p2, mode="r")
        ]
        collated_mags = collate_feature_data_mags(mags)

        # compare directories
        expected = self.get_data_path("collated/feature_data")
        dircmp = filecmp.dircmp(collated_mags.path, expected)
        self.assertListEqual(
            [
                "24dee6fe-9b84-45bb-8145-de7b092533a1.fasta",
                "fb0bc871-04f6-486b-a10e-8e0cb66f8de3.fasta"
            ],
            dircmp.common
        )
