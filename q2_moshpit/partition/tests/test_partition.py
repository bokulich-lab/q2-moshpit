# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from unittest.mock import patch
from q2_moshpit.partition.partition import partition_mags, collate_mags
from qiime2.plugin.testing import TestPluginBase
from q2_types_genomics.per_sample_data._format import (
    MultiMAGSequencesDirFmt, MultiFASTADirectoryFormat
)
from q2_types_genomics.feature_data._format import (
    MAGSequencesDirFmt
)


class TestHelpers(TestPluginBase):
    package = 'q2_moshpit.partition.tests'

    @patch('qiime2.util.duplicate')
    @patch('q2_moshpit.partition.utils._validate_num_partitions')
    @patch('q2_moshpit.partition.utils._partition_sample_data_mags')
    @patch('q2_moshpit.partition.utils._partition_feature_data_mags')
    def test_partition_mags_sample(
        self,
        mock_feature,
        mock_sample,
        mock_val_num,
        mock_duplicate
    ):
        mags = MultiMAGSequencesDirFmt()
        partition_mags(mags)
        mock_sample.assert_called_once_with(mags, None)
        mock_feature.assert_not_called()

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

    def test_partition_mags_neither(self):
        mags = [MultiFASTADirectoryFormat()]
        with self.assertRaisesRegex(
            ValueError,
            "--i-mags is neither"
        ):
            partition_mags(mags)

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

    def test_collate_mags_neither(self):
        mags = [MultiFASTADirectoryFormat()]
        with self.assertRaisesRegex(
            ValueError,
            "--i-mags is neither"
        ):
            collate_mags(mags)
