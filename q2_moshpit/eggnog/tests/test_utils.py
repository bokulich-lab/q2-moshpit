# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
import tempfile
from unittest.mock import patch, call
from qiime2.plugin.testing import TestPluginBase
from q2_types.reference_db import HmmerDirFmt
from .._dbs import _validate_taxon_id, _try_wget, _download_and_build_hmm_db


class TestEggnogUtils(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    def test_validate_taxon_id_invalid(self):
        path_to_data = self.get_data_path('build_eggnog_diamond_db/')
        with self.assertRaisesRegex(ValueError, "'0' is not valid taxon ID. "):
            _validate_taxon_id(path_to_data, 0)

    def test_validate_taxon_id_valid(self):
        path_to_data = self.get_data_path('build_eggnog_diamond_db/')
        _validate_taxon_id(path_to_data, 2)

    @patch("subprocess.run")
    def test_try_wget(self, mock_run):
        _try_wget("foo.txt", "www.fake_url.com", "Download error.")
        mock_run.assert_called_once_with(
            ["wget", "-O", "foo.txt", "www.fake_url.com"], check=True
        )

    @patch("subprocess.run")
    def test_try_wget_exception(self, mock_run):
        mock_run.side_effect = subprocess.CalledProcessError(
            1, ["wget", "-O", "foo.txt", "www.fake_url.com"]
        )
        with self.assertRaisesRegex(Exception, "Download error: 1"):
            _try_wget("foo.txt", "www.fake_url.com", "Download error")

    @patch("q2_moshpit.eggnog._utils._try_wget")
    @patch("subprocess.run")
    @patch.object(tempfile, 'TemporaryDirectory')
    def test_download_and_build_hmm_db(self, tmpdir, mock_run, mock_wet):
        tmp = self.get_data_path('hmmer')
        taxon_id = 1
        tmpdir.return_value.__enter__.return_value = tmp

        hmmer_db = _download_and_build_hmm_db(taxon_id)

        self.assertIsInstance(hmmer_db, HmmerDirFmt)
        mock_wet.assert_called_once_with(
            f"{tmp}/{taxon_id}_hmms.tar.gz",
            "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/"
            f"{taxon_id}/{taxon_id}_hmms.tar.gz",
            "Error during HMMER database download"
        )
        mock_run.assert_has_calls([
            call(
                ["tar", "zxf", f"{taxon_id}_hmms.tar.gz"],
                check=True,
                cwd=tmp
            ),
            call(
                ["hmmpress", f"{str(hmmer_db)}/{taxon_id}.hmm"],
                check=True
            )
        ])
