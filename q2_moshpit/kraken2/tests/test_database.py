# ----------------------------------------------------------------------------
# Copyright (c) 2023-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import io
import os
import shutil
import tarfile
import tempfile
import unittest
from copy import deepcopy
from requests.exceptions import ConnectionError
from subprocess import CalledProcessError
from tempfile import TemporaryDirectory
from unittest.mock import patch, ANY, call, Mock, MagicMock

from q2_types.feature_data import DNAFASTAFormat
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import moshpit

from q2_moshpit.kraken2.database import (
    _fetch_taxonomy, _fetch_libraries,
    _add_seqs_to_library, _build_kraken2_database, _move_db_files, _build_bracken_database, _find_latest_db,
    _fetch_db_collection, S3_COLLECTIONS_URL, _build_dbs_from_seqs, _fetch_prebuilt_dbs
)
from q2_types_genomics.kraken2 import Kraken2DBDirectoryFormat, BrackenDBDirectoryFormat


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestKraken2Database(TestPluginBase):
    package = "q2_moshpit.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.kraken2_db_dir = 'fake/db/dir'
        self.kwargs = {
            'threads': 2, 'fast_build': True,
            'kmer_len': 31, 'use_ftp': False,
            'max_db_size': 1000, 'load_factor': 0.5
        }
        self.s3_response = b'''
            <ListBucketResult>
                <Contents>
                    <Key>kraken/k2_viral_20201202.tar.gz</Key>
                    <LastModified>2020-12-09T01:38:22.000Z</LastModified>
                </Contents>
                <Contents>
                    <Key>kraken/k2_viral_20230314.tar.gz</Key>
                    <LastModified>2023-03-22T01:29:11.000Z</LastModified>
                </Contents>
            </ListBucketResult>
        '''

        self.temp_dir = tempfile.mkdtemp()
        self.temp_tar = os.path.join(self.temp_dir, 'temp.tar.gz')

        with tarfile.open(self.temp_tar, "w:gz") as tar:
            data = io.BytesIO(b"sample data")
            tarinfo = tarfile.TarInfo(name="sample.txt")
            tarinfo.size = len(data.getbuffer())
            tar.addfile(tarinfo, data)

        with open(self.temp_tar, "rb") as f:
            self.tar_chunks = [
                chunk for chunk in iter(lambda: f.read(8192), b"")
            ]

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_fetch_taxonomy(self, p1):
        _fetch_taxonomy(self.kraken2_db_dir, threads=3, use_ftp=True)

        exp_cmd = [
            "kraken2-build", "--download-taxonomy",
            "--threads", "3", "--db", self.kraken2_db_dir, "--use-ftp"
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch(
        'q2_moshpit.kraken2.database.run_command',
        side_effect=CalledProcessError(123, 'cmd')
    )
    def test_fetch_taxonomy_exception(self, p1):
        with self.assertRaisesRegex(
                Exception,
                "An error was encountered .* downloading taxonomy, "
                r"\(return code 123\), please inspect .*"
        ):
            _fetch_taxonomy(self.kraken2_db_dir, threads=3, use_ftp=True)

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_fetch_libraries_skip(self, p1):
        all_kwargs = deepcopy(self.kwargs)
        all_kwargs['library_exists'] = 'skip'
        libraries = ['plasmid', 'human']

        with TemporaryDirectory() as tmp_dir:
            for lib in libraries:
                os.makedirs(os.path.join(tmp_dir, 'library', lib))
            open(os.path.join(
                tmp_dir, 'library', libraries[0], 'library.fna'), 'w'
            ).close()

            _fetch_libraries(
                tmp_dir, libraries=libraries, all_kwargs=all_kwargs
            )

        exp_cmd = [
            "kraken2-build", "--download-library", libraries[1],
            "--threads", "2", "--db", tmp_dir
        ]
        p1.assert_called_once_with(
            cmd=exp_cmd, verbose=True
        )

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_fetch_libraries_refetch(self, p1):
        all_kwargs = deepcopy(self.kwargs)
        all_kwargs['library_exists'] = 'refetch'
        libraries = ['plasmid', 'human']

        _fetch_libraries(
            self.kraken2_db_dir, libraries=libraries, all_kwargs=all_kwargs
        )

        base_cmd = ["kraken2-build", "--download-library"]
        exp_common_args = ["--threads", "2", "--db", self.kraken2_db_dir]
        exp_cmds = [
            [*base_cmd, libraries[0], *exp_common_args],
            [*base_cmd, libraries[1], *exp_common_args]
        ]
        p1.assert_has_calls([
            call(cmd=exp_cmds[0], verbose=True),
            call(cmd=exp_cmds[1], verbose=True)
        ])

    @patch(
        'q2_moshpit.kraken2.database.run_command',
        side_effect=CalledProcessError(123, 'cmd')
    )
    def test_fetch_libraries_exception(self, p1):
        with self.assertRaisesRegex(
                Exception,
                "An error was encountered .* downloading the 'human' "
                r"library, \(return code 123\), please inspect .*"
        ):

            _fetch_libraries(
                self.kraken2_db_dir, libraries=['human'], all_kwargs=self.kwargs
            )

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_add_seqs_to_library(self, p1):
        seqs = DNAFASTAFormat(self.get_data_path('mags/samp1/bin1.fa'), 'r')

        _add_seqs_to_library(self.kraken2_db_dir, seqs=seqs, no_masking=True)

        exp_cmd = [
            "kraken2-build", "--add-to-library", str(seqs.path),
            "--db", self.kraken2_db_dir, "--no-mask"
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch(
        'q2_moshpit.kraken2.database.run_command',
        side_effect=CalledProcessError(123, 'cmd')
    )
    def test_add_seqs_to_library_exception(self, p1):
        seqs = DNAFASTAFormat(self.get_data_path('mags/samp1/bin1.fa'), 'r')

        with self.assertRaisesRegex(
                Exception,
                "An error was encountered .* adding sequences to the "
                r"library, \(return code 123\), please inspect .*"
        ):
            _add_seqs_to_library(self.kraken2_db_dir, seqs=seqs, no_masking=True)

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_build_kraken2_database(self, p1):
        _build_kraken2_database(self.kraken2_db_dir, all_kwargs=self.kwargs)

        exp_cmd = [
            "kraken2-build", "--build", "--db", self.kraken2_db_dir,
            "--threads", "2", "--fast-build", "--kmer-len", "31",
            "--load-factor", "0.5", "--max-db-size", "1000"
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_build_kraken2_database_no_max_db(self, p1):
        all_kwargs = deepcopy(self.kwargs)
        all_kwargs['max_db_size'] = 0

        _build_kraken2_database(self.kraken2_db_dir, all_kwargs=all_kwargs)

        exp_cmd = [
            "kraken2-build", "--build", "--db", self.kraken2_db_dir,
            "--threads", "2", "--fast-build", "--kmer-len", "31",
            "--load-factor", "0.5"
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch(
        'q2_moshpit.kraken2.database.run_command',
        side_effect=CalledProcessError(123, 'cmd')
    )
    def test_build_kraken2_database_exception(self, p1):
        with self.assertRaisesRegex(
                Exception,
                "An error was encountered .* building the database, "
                r"\(return code 123\), please inspect .*"
        ):
            _build_kraken2_database(
                self.kraken2_db_dir, all_kwargs=self.kwargs
            )

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_build_bracken_database(self, p1):
        _build_bracken_database(
            kraken2_db_dir=self.kraken2_db_dir, threads=2,
            kmer_len=31, read_len=150
        )

        exp_cmd = [
            "bracken-build", "-d", self.kraken2_db_dir,
            "-t", "2", "-k", "31", "-l", "150"
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch(
        'q2_moshpit.kraken2.database.run_command',
        side_effect=CalledProcessError(123, 'cmd')
    )
    def test_build_bracken_database_exception(self, p1):
        with self.assertRaisesRegex(
                Exception,
                "An error was encountered while building the Bracken "
                r"database, \(return code 123\), please inspect .+"
        ):
            _build_bracken_database(
                kraken2_db_dir=self.kraken2_db_dir, threads=2,
                kmer_len=31, read_len=150
            )

    def test_find_latest_db(self):
        response = Mock(content=self.s3_response)

        obs_db = _find_latest_db('viral', response)
        exp_db = 'kraken/k2_viral_20230314.tar.gz'
        self.assertEqual(obs_db, exp_db)

    def test_find_latest_db_empty(self):
        response = Mock(content=b'''<ListBucketResult></ListBucketResult>''')

        with self.assertRaisesRegex(
            ValueError, r'No databases were found.+'
        ):
            _find_latest_db('viral', response)

    @patch("requests.get")
    @patch("tarfile.open")
    @patch(
        "q2_moshpit.kraken2.database._find_latest_db",
        return_value="kraken/k2_viral.tar.gz"
    )
    def test_fetch_db_collection_success(
            self, mock_find, mock_tarfile_open, mock_requests_get
    ):
        mock_requests_get.side_effect = [
            MagicMock(status_code=200),
            MagicMock(
                status_code=200,
                iter_content=lambda chunk_size: self.tar_chunks
            )
        ]
        mock_tarfile_open.return_value.__enter__.return_value = MagicMock()

        _fetch_db_collection("viral", "/tmp")

        mock_requests_get.has_calls([
            call(S3_COLLECTIONS_URL),
            call(f"{S3_COLLECTIONS_URL}/kraken/k2_viral.tar.gz", stream=True)
        ])
        mock_tarfile_open.assert_called_once_with(
            "/tmp/k2_viral.tar.gz", "r:gz"
        )
        mock_find.assert_called_once_with("viral", ANY)

    @patch('requests.get')
    def test_fetch_db_collection_connection_error(self, mock_get):
        mock_get.side_effect = ConnectionError("Some error.")
        with self.assertRaisesRegex(
                ValueError, r".+The error was\: Some error\."
        ):
            _fetch_db_collection('my_collection', '/tmp')

    @patch('requests.get')
    def test_fetch_db_collection_status_non200(self, mock_get):
        mock_get.return_value = Mock(status_code=404)
        with self.assertRaisesRegex(
                ValueError, r".+Status code was\: 404"
        ):
            _fetch_db_collection('my_collection', '/tmp')

    def test_move_db_files(self):
        with TemporaryDirectory() as tmp_dir:
            fake_src = os.path.join(tmp_dir, 'fake_src')
            fake_dest = os.path.join(tmp_dir, 'fake_dest')
            fake_files = ['fake_db.k2d', 'fake_file.k2d', 'other.file']

            os.makedirs(fake_src)
            os.makedirs(fake_dest)

            for f in fake_files:
                open(os.path.join(tmp_dir, 'fake_src', f), 'w').close()

            _move_db_files(fake_src, fake_dest)

            for f in fake_files[:2]:
                self.assertTrue(os.path.exists(os.path.join(fake_dest, f)))

    @patch("q2_moshpit.kraken2.database._fetch_taxonomy")
    @patch("q2_moshpit.kraken2.database._add_seqs_to_library")
    @patch("q2_moshpit.kraken2.database._build_kraken2_database")
    @patch("q2_moshpit.kraken2.database._build_bracken_database")
    @patch("q2_moshpit.kraken2.database._move_db_files")
    def test_build_dbs_from_seqs(
            self, mock_move, mock_bracken, mock_kraken,
            mock_add_seqs, mock_fetch_tax
    ):
        bracken_db, kraken2_db = MagicMock(), MagicMock()
        seqs, tmp_dir = ["seq1", "seq2"], "/tmp"
        common_args = {
            "threads": 1, "use_ftp": False, "no_masking": False,
            "read_len": [100, 150], "kmer_len": 35
        }

        _build_dbs_from_seqs(
            bracken_db, kraken2_db, seqs, tmp_dir, common_args
        )

        mock_fetch_tax.assert_called_once_with(
            db_dir=tmp_dir, threads=1, use_ftp=False
        )
        mock_add_seqs.assert_has_calls([
            call(db_dir=tmp_dir, seqs="seq1", no_masking=False),
            call(db_dir=tmp_dir, seqs="seq2", no_masking=False)
        ])
        mock_kraken.assert_called_once_with(
            db_dir=tmp_dir, all_kwargs=common_args
        )
        mock_bracken.assert_has_calls([
            call(kraken2_db_dir=tmp_dir, threads=1,
                 kmer_len=35, read_len=100),
            call(kraken2_db_dir=tmp_dir, threads=1,
                 kmer_len=35, read_len=150)
        ])
        mock_move.has_calls([
            call(tmp_dir, str(kraken2_db.path), extension="k2d"),
            call(tmp_dir, str(bracken_db.path), extension="kmer_distrib")
        ])

    @patch("q2_moshpit.kraken2.database._fetch_db_collection")
    @patch("q2_moshpit.kraken2.database._move_db_files")
    def test_fetch_prebuilt_dbs(self, mock_move, mock_fetch):
        bracken_db = MagicMock(path="/path/to/bracken_db")
        kraken2_db = MagicMock(path="/path/to/kraken2_db")

        _fetch_prebuilt_dbs(bracken_db, kraken2_db, "some_collection", "/tmp")

        mock_fetch.assert_called_once_with(
            collection="some_collection", tmp_dir="/tmp"
        )
        mock_move.assert_has_calls([
            call("/tmp", str(kraken2_db.path), extension="k2d"),
            call("/tmp", str(bracken_db.path), extension="kmer_distrib")
        ])

    @patch("tempfile.TemporaryDirectory", return_value=MockTempDir())
    @patch("q2_moshpit.kraken2.database.Kraken2DBDirectoryFormat")
    @patch("q2_moshpit.kraken2.database.BrackenDBDirectoryFormat")
    @patch("q2_moshpit.kraken2.database._fetch_prebuilt_dbs")
    def test_build_kraken_db_action_with_prebuilt(
            self, mock_fetch, mock_bracken, mock_kraken, mock_tmp
    ):
        fake_kraken_dir_fmt = Kraken2DBDirectoryFormat(
            self.get_data_path('db'), 'r'
        )
        mock_kraken.return_value = fake_kraken_dir_fmt
        fake_bracken_dir_fmt = BrackenDBDirectoryFormat(
            self.get_data_path('bracken-db'), 'r'
        )
        mock_bracken.return_value = fake_bracken_dir_fmt

        moshpit.actions.build_kraken_db(collection="viral")

        mock_fetch.assert_called_once_with(
            fake_bracken_dir_fmt, fake_kraken_dir_fmt,
            "viral", str(mock_tmp.return_value.name)
        )

    @patch("tempfile.TemporaryDirectory", return_value=MockTempDir())
    @patch("q2_moshpit.kraken2.database.Kraken2DBDirectoryFormat")
    @patch("q2_moshpit.kraken2.database.BrackenDBDirectoryFormat")
    @patch("q2_moshpit.kraken2.database._build_dbs_from_seqs")
    def test_build_kraken_db_action_with_seqs(
            self, mock_build, mock_bracken, mock_kraken, mock_tmp
    ):
        seqs = Artifact.import_data(
            'FeatureData[Sequence]', self.get_data_path("seqs")
        )
        fake_kraken_dir_fmt = Kraken2DBDirectoryFormat(
            self.get_data_path('db'), 'r'
        )
        mock_kraken.return_value = fake_kraken_dir_fmt
        fake_bracken_dir_fmt = BrackenDBDirectoryFormat(
            self.get_data_path('bracken-db'), 'r'
        )
        mock_bracken.return_value = fake_bracken_dir_fmt

        moshpit.actions.build_kraken_db(
            seqs=[seqs], threads=2, fast_build=True
        )

        exp_common_args = {
            'threads': 2, 'kmer_len': 35, 'minimizer_len': 31,
            'minimizer_spaces': 7, 'no_masking': False, 'max_db_size': 0,
            'use_ftp': False, 'load_factor': 0.7, 'fast_build': True,
            'read_len': [50, 75, 100, 150, 200, 250, 300],
            'kraken2_db': fake_kraken_dir_fmt,
            'bracken_db': fake_bracken_dir_fmt,
            'tmp': str(mock_tmp.return_value.name)
        }
        mock_build.assert_called_once_with(
            fake_bracken_dir_fmt, fake_kraken_dir_fmt,
            [ANY], str(mock_tmp.return_value.name),
            exp_common_args
        )

    @patch("tempfile.TemporaryDirectory", return_value=MockTempDir())
    def test_build_kraken_db_action_with_error(self, mock_tmp):
        with self.assertRaisesRegex(
            ValueError, r"You need to either provide a list .+"
        ):
            moshpit.actions.build_kraken_db()


if __name__ == "__main__":
    unittest.main()
