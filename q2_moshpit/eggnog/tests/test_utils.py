# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
from unittest.mock import patch, Mock
from qiime2.plugin.testing import TestPluginBase
from .._dbs import _validate_taxon_id, _try_wget


class TestFetchDB(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    def test_validate_taxon_id_invalid(self):
        # Init input data
        path_to_data = self.get_data_path('build_eggnog_diamond_db/')

        # Call function exception error since taxon 0 is invalid
        with self.assertRaisesRegex(
            ValueError,
            "'0' is not valid taxon ID. "
        ):
            _validate_taxon_id(path_to_data, 0)

    def test_validate_taxon_id_valid(self):
        # Init input data
        path_to_data = self.get_data_path('build_eggnog_diamond_db/')
        _validate_taxon_id(path_to_data, 2)

    @patch("subprocess.run")
    def test_try_wget(self, mock_run):
        _try_wget("foo.txt", "www.fake_url.com", "Download error.")
        mock_run.assert_called_once_with(
            ["wget", "-O", "foo.txt", "www.fake_url.com"],
            check=True
        )

    @patch("subprocess.run")
    def test_try_wget_exception(self, mock_run):
        mock_run.side_effect = subprocess.CalledProcessError(
            1,
            ["wget", "-O", "foo.txt", "www.fake_url.com"]
        )
        with self.assertRaisesRegex(
            Exception,
            "Download error: 1"
        ):
            _try_wget("foo.txt", "www.fake_url.com", "Download error")
