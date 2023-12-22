# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from unittest.mock import patch
from qiime2.plugin.testing import TestPluginBase
from .._method import fetch_eggnog_db


class TestFetchDB(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    @patch("subprocess.run")
    def test_fetch_eggnog_db(self, subp_run):
        # Call function. Patching will make sure nothing is
        # actually ran
        eggnog_db = fetch_eggnog_db()

        # Check that command was called in the expected way
        cmd = [
            "download_eggnog_data.py", "-y", "-D",
            "--data_dir", str(eggnog_db)
        ]
        subp_run.assert_called_once_with(cmd, check=True)
