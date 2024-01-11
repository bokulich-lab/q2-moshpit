# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
from unittest.mock import patch, call
from qiime2.plugin.testing import TestPluginBase
from .._dbs import fetch_eggnog_db, fetch_eggnog_proteins, fetch_diamond_db


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

    @patch("subprocess.run")
    def test_fetch_diamond_db(self, subp_run):
        # Call function. Patching will make sure nothing is
        # actually ran
        diamond_db = fetch_diamond_db()
        path_out = os.path.join(str(diamond_db), "ref_db.dmnd.gz")

        # Check that command was called in the expected way
        first_call = call(
            [
                "wget", "-e", "robots=off", "-O", f"{path_out}",
                "http://eggnogdb.embl.de/download/emapperdb-5.0.2/"
                "eggnog_proteins.dmnd.gz"
            ],
            check=True
        )
        second_call = call(
            ["gunzip", f"{path_out}"],
            check=True,
        )

        # Check that commands are ran as expected
        subp_run.assert_has_calls([first_call, second_call], any_order=False)

    @patch("subprocess.run")
    def test_fetch_eggnog_fasta(self, subp_run):
        # Call function. Patching will make sure nothing is
        # actually ran
        eggnog_fa = fetch_eggnog_proteins()
        fasta_file = os.path.join(str(eggnog_fa), "e5.proteomes.faa")
        taxonomy_file = os.path.join(str(eggnog_fa), "e5.taxid_info.tsv")

        # Check that command was called in the expected way
        first_call = call(
            [
                "wget", "-e", "robots=off", "-O", f"{fasta_file}",
                "http://eggnog5.embl.de/download/eggnog_5.0/e5.proteomes.faa"
            ],
            check=True
        )
        second_call = call(
            [
                "wget", "-e", "robots=off", "-O", f"{taxonomy_file}",
                "http://eggnog5.embl.de/download/eggnog_5.0/e5.taxid_info.tsv"
            ],
            check=True,
        )

        # Check that commands are ran as expected
        subp_run.assert_has_calls([first_call, second_call], any_order=False)
