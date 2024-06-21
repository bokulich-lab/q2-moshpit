# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
from unittest.mock import patch, call
from qiime2.plugin.testing import TestPluginBase
from .._utils import (
    _try_wget, _download_and_build_hmm_db, _download_fastas_into_hmmer_db,
    _merge_hmms_and_write_idmap
)
from q2_moshpit.eggnog._format import EggnogHmmerIdmapDirectoryFmt
from q2_types.genome_data import ProteinsDirectoryFormat
from q2_types.profile_hmms import (
    PressedProfileHmmsDirectoryFmt, ProteinMultipleProfileHmmDirectoryFmt
)


class TestEggnogUtils(TestPluginBase):
    package = "q2_moshpit.eggnog.tests"

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
    @patch("tempfile.TemporaryDirectory")
    def test_download_and_build_hmm_db(self, tmpdir, mock_run, mock_wet):
        tmp = self.get_data_path('hmmer/hmms')
        taxon_id = 1
        tmpdir.return_value.__enter__.return_value = tmp

        idmap, hmm_db, pressed_hmm_db = _download_and_build_hmm_db(taxon_id)

        self.assertIsInstance(idmap, EggnogHmmerIdmapDirectoryFmt)
        self.assertIsInstance(hmm_db, ProteinMultipleProfileHmmDirectoryFmt)
        self.assertIsInstance(pressed_hmm_db, PressedProfileHmmsDirectoryFmt)

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
                ["hmmpress", f"{str(pressed_hmm_db)}/{taxon_id}.hmm"],
                check=True
            )
        ])

    @patch("q2_moshpit.eggnog._utils._try_wget")
    @patch("subprocess.run")
    @patch("tempfile.TemporaryDirectory")
    def test_download_fastas_into_hmmer_db(self, tmpdir, mock_run, mock_wet):
        tmp = self.get_data_path("hmmer/fastas")
        taxon_id = 1
        tmpdir.return_value.__enter__.return_value = tmp
        fastas = _download_fastas_into_hmmer_db(1)

        self.assertIsInstance(fastas, ProteinsDirectoryFormat)

        mock_wet.assert_called_once_with(
            f"{tmp}/{taxon_id}_raw_algs.tar",
            "http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level/"
            f"{taxon_id}/{taxon_id}_raw_algs.tar",
            "Error downloading FASTA files"
        )
        mock_run.assert_called_once_with(
            ["tar", "xf", f"{taxon_id}_raw_algs.tar"],
            check=True,
            cwd=tmp
        )

    def test_merge_hmms_and_write_idmap(self):
        hmms = self.get_data_path("hmmer/hmms")
        taxon_id = 1

        merged_hmms_obj = ProteinMultipleProfileHmmDirectoryFmt()
        idmap_obj = EggnogHmmerIdmapDirectoryFmt()

        hmms_merged_p = f"{str(merged_hmms_obj)}/{taxon_id}.hmm"
        idmap_p = f"{str(idmap_obj)}/{taxon_id}.hmm.idmap"

        _merge_hmms_and_write_idmap(hmms_merged_p, idmap_p, taxon_id, hmms)

        merged_hmms_obj.validate()
        idmap_obj.validate()
