# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pdb
import filecmp
from unittest.mock import patch
from q2_moshpit.partition.partition import (
    partition_sample_data_mags, partition_feature_data_mags
)
from qiime2.plugin.testing import TestPluginBase
from q2_types_genomics.per_sample_data._format import (
    MultiMAGSequencesDirFmt, MultiFASTADirectoryFormat
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
        p = self.get_data_path("sample_data")
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

    @patch('qiime2.util.duplicate')
    @patch('q2_moshpit.partition.utils._validate_num_partitions')
    @patch('q2_moshpit.partition.utils._partition_sample_data_mags')
    @patch('q2_moshpit.partition.utils._partition_feature_data_mags')
    def test_partition_mags_feature(
        self,
        mock_feature,
        mock_sample,
        mock_val_num,
        mock_duplicate
    ):
        mags = MAGSequencesDirFmt()
        partition_mags(mags)
        mock_feature.assert_called_once_with(mags, None)
        mock_sample.assert_not_called()

    @patch('q2_moshpit.partition.utils._collate_sample_data_mags')
    @patch('q2_moshpit.partition.utils._collate_feature_data_mags')
    def test_collate_mags_sample(self, mock_feature, mock_sample):
        mags = [MultiMAGSequencesDirFmt()]
        collate_mags(mags)
        mock_feature.assert_not_called()
        mock_sample.assert_called_once_with(mags)

    @patch('q2_moshpit.partition.utils._collate_sample_data_mags')
    @patch('q2_moshpit.partition.utils._collate_feature_data_mags')
    def test_collate_mags_features(self, mock_feature, mock_sample):
        mags = [MAGSequencesDirFmt()]
        collate_mags(mags)
        mock_sample.assert_not_called()
        mock_feature.assert_called_once_with(mags)
