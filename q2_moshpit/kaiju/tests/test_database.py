# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
import unittest
from unittest.mock import patch, Mock

import pandas as pd
from bs4 import BeautifulSoup
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.kaiju.database import (
    _fetch_and_extract_db, _find_latest_db_url, _find_all_dbs,
    fetch_kaiju_db, CHUNK_SIZE, ERR_MSG, KAIJU_SERVER_URL
)
from requests.exceptions import ConnectionError, RequestException

from q2_types_genomics.kaiju import KaijuDBDirectoryFormat


class TestDatabaseFunctions(TestPluginBase):
    package = 'q2_moshpit.kaiju.tests'

    @patch("requests.get")
    @patch("q2_moshpit.kaiju.database.tqdm")
    @patch("tarfile.open")
    @patch("os.remove")
    def test_fetch_and_extract_db(
            self, mock_remove, mock_tarfile_open,
            mock_progress, mock_requests
    ):
        response = mock_requests.return_value
        response.headers = {"content-length": 1024}
        response.iter_content.return_value = [b"test"] * 1024
        mock_tar = Mock()
        mock_tarfile_open.return_value.__enter__.return_value = mock_tar

        with tempfile.TemporaryDirectory() as tmpdir:
            _fetch_and_extract_db("http://a/b/db.tar.gz", tmpdir)
            db_path = os.path.join(tmpdir, "db.tar.gz")

            mock_progress.assert_called_with(
                desc='Downloading the "db.tar.gz" database',
                total=1024,
                unit="B",
                unit_scale=True,
                unit_divisor=1024
            )
            response.iter_content.assert_called_with(chunk_size=CHUNK_SIZE)
            mock_tarfile_open.assert_called_with(db_path, "r:gz")
            mock_tar.extractall.assert_called_with(path=tmpdir)
            mock_remove.assert_called_with(db_path)
            mock_requests.assert_called_with(
                "http://a/b/db.tar.gz", stream=True
            )

    @patch("requests.get", side_effect=ConnectionError("some error"))
    def test_fetch_and_extract_db_exception(
            self, mock_requests
    ):
        exp_error = ERR_MSG.format("some error")
        with self.assertRaisesRegex(Exception, exp_error):
            with tempfile.TemporaryDirectory() as tmpdir:
                _fetch_and_extract_db("http://a/b/db.tar.gz", tmpdir)

                mock_requests.assert_called_with(
                    "http://a/b/db.tar.gz", stream=True
                )

    def test_find_latest_db_url(self):
        databases = [
            ('nr_euk 2021-02-24 (61GB)',
             'https://hello.com/nr_euk_2021-02-24.tar.gz'),
            ('nr 2021-02-26 (52GB)',
             'https://hello.com/nr_2021-02-26.tar.gz'),
            ('nr_euk 2022-01-11 (60GB)',
             'https://hello.com/nr_euk_2022-01-11.tar.gz')
        ]
        sidebox_element = BeautifulSoup(
            '<html><body>{}</body></html>'.format(
                ''.join('<a href={}>{}</a>'.format(d[1], d[0])
                        for d in databases)
            ), 'html.parser')
        url = _find_latest_db_url(
            database_type='nr_euk',
            sidebox_element=sidebox_element,
            url='https://test.com'
        )
        self.assertEqual(url, 'https://hello.com/nr_euk_2022-01-11.tar.gz')

    def test_find_all_dbs(self):
        databases = ['nr_euk 2021-02-24 (61GB)', 'nr 2021-02-26 (52GB)']
        sidebox_element = BeautifulSoup(
            '<html><body>{}</body></html>'.format(
                ''.join('<a>{}</a>'.format(d) for d in databases)
            ), 'html.parser')
        df = _find_all_dbs(sidebox_element)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertListEqual(
            df.index.tolist(),
            ['nr_euk 2021-02-24 (61GB)', 'nr 2021-02-26 (52GB)']
        )
        self.assertListEqual(
            df['Date'].tolist(),
            [pd.to_datetime('2021-02-24'), pd.to_datetime('2021-02-26')]
        )

    @patch("requests.get")
    @patch("q2_moshpit.kaiju.database._fetch_and_extract_db")
    def test_fetch_kaiju_db(self, mock_fetch, mock_requests):
        databases = [
            ('nr_euk 2021-02-24 (61GB)',
             'https://hello.com/nr_euk_2021-02-24.tar.gz'),
            ('nr 2021-02-26 (52GB)',
             'https://hello.com/nr_2021-02-26.tar.gz'),
            ('nr_euk 2022-01-11 (60GB)',
             'https://hello.com/nr_euk_2022-01-11.tar.gz')
        ]
        mock_requests.return_value = Mock(
            content='<html><body><div id="sidebox_db">{}</div></body></html>'
            .format(
                ''.join('<a href={}>{}</a>'.format(d[1], d[0])
                        for d in databases)
            )
        )

        obs_db = fetch_kaiju_db('nr_euk')
        self.assertIsInstance(obs_db, KaijuDBDirectoryFormat)
        mock_requests.assert_called_with(KAIJU_SERVER_URL)
        mock_fetch.assert_called_with(
            'https://hello.com/nr_euk_2022-01-11.tar.gz',
            str(obs_db.path)
        )

    @patch("requests.get", side_effect=RequestException("some error"))
    def test_fetch_kaiju_db_exception(self, mock_requests):
        with self.assertRaisesRegex(
                Exception, ERR_MSG.format("some error")
        ):
            fetch_kaiju_db('nr_euk')

        mock_requests.assert_called_with(KAIJU_SERVER_URL)


if __name__ == "__main__":
    unittest.main()
