# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import tempfile
import pandas as pd
from q2_moshpit.busco.busco import _run_busco, _busco_helper
from unittest.mock import patch, ANY, call
from qiime2.plugin.testing import TestPluginBase
from q2_types.per_sample_sequences._format import MultiMAGSequencesDirFmt


class TestBUSCO(TestPluginBase):
    package = "q2_moshpit.busco.tests"

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.mags = MultiMAGSequencesDirFmt(
            path=self.get_data_path('mags'),
            mode="r",
        )

    def _prepare_summaries(self):
        for s in ['1', '2']:
            os.makedirs(os.path.join(self.temp_dir.name, f"sample{s}"))
            shutil.copy(
                self.get_data_path(f'summaries/batch_summary_{s}.txt'),
                os.path.join(
                    self.temp_dir.name, f"sample{s}", 'batch_summary.txt'
                )
            )

    @patch('q2_moshpit.busco.busco.run_command')
    def test_run_busco(self, mock_run):
        self._prepare_summaries()

        obs = _run_busco(
            output_dir=self.temp_dir.name,
            mags=self.mags,
            params=['--lineage_dataset', 'bacteria_odb10', '--cpu', '7']
        )
        exp = {
            'sample1': f"{self.temp_dir.name}/sample1/batch_summary.txt",
            'sample2': f"{self.temp_dir.name}/sample2/batch_summary.txt",
        }

        self.assertDictEqual(obs, exp)
        mock_run.assert_has_calls([
            call(
                ['busco', '--lineage_dataset', 'bacteria_odb10',
                 '--cpu', '7', '--in', self.get_data_path('mags/sample1'),
                 '--out_path', self.temp_dir.name, '-o', 'sample1'],
            ),
            call(
                ['busco', '--lineage_dataset', 'bacteria_odb10',
                 '--cpu', '7', '--in', self.get_data_path('mags/sample2'),
                 '--out_path', self.temp_dir.name, '-o', 'sample2'],
            )
        ])

    @patch('q2_moshpit.busco.busco._run_busco')
    @patch('q2_moshpit.busco.busco._get_mag_lengths')
    def test_busco_helper(self, mock_len, mock_run):
        self._prepare_summaries()
        mock_run.return_value = {
            'sample1': f"{self.temp_dir.name}/sample1/batch_summary.txt",
            'sample2': f"{self.temp_dir.name}/sample2/batch_summary.txt",
        }
        mock_len.return_value = pd.Series(
            {
                'ab23d75d-547d-455a-8b51-16b46ddf7496': 3177889,
                '0e514d88-16c4-4273-a1df-1a360eb2c823': 19625518,
                '8098e3a8-df4a-46af-83e2-6c2443d74cb9': 912062,
                'd4637408-80ab-49d5-ab32-c66509c3a544': 2678999,
                '503c2f56-3e4f-4ce7-9b61-b63bc7fe0592': 21035714
            }, name='length'
        )

        obs = _busco_helper(
            self.mags, ['--lineage_dataset', 'bacteria_odb10']
        )
        exp = pd.read_csv(
            self.get_data_path('summaries/all_renamed_with_lengths.csv'),
        )

        pd.testing.assert_frame_equal(obs, exp)
        mock_run.assert_called_with(
            output_dir=ANY, mags=self.mags,
            params=['--lineage_dataset', 'bacteria_odb10']
        )
