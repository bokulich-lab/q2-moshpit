# ----------------------------------------------------------------------------
# Copyright (c) 2023-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
import unittest
from copy import deepcopy
from subprocess import CalledProcessError
from tempfile import TemporaryDirectory
from unittest.mock import patch, ANY, call

from q2_types.feature_data import DNAFASTAFormat
from qiime2 import Artifact
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugins import moshpit

from q2_moshpit.kraken2.database import (
    _build_standard_db, _fetch_taxonomy, _fetch_libraries,
    _add_seqs_to_library, _build_database, _move_db_files
)
from q2_types_genomics.kraken2 import Kraken2DBDirectoryFormat


class MockTempDir(tempfile.TemporaryDirectory):
    pass


class TestKraken2Database(TestPluginBase):
    package = "q2_moshpit.kraken2.tests"

    def setUp(self):
        super().setUp()
        self.db_dir = 'fake/db/dir'
        self.kwargs = {
            'threads': 2, 'fast_build': True,
            'kmer_len': 31, 'use_ftp': False,
            'max_db_size': 1000, 'load_factor': 0.5
        }

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_build_standard_db(self, p1):
        _build_standard_db(self.db_dir, self.kwargs)

        exp_cmd = [
            "kraken2-build", "--standard", "--db", self.db_dir,
            "--threads", "2", "--fast-build", "--max-db-size", "1000"
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch(
        'q2_moshpit.kraken2.database.run_command',
        side_effect=CalledProcessError(123, 'cmd')
    )
    def test_build_standard_db_exception(self, p1):
        with self.assertRaisesRegex(
                Exception,
                "An error was encountered .* standard library, "
                r"\(return code 123\), please inspect .*"
        ):
            _build_standard_db(self.db_dir, self.kwargs)

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_fetch_taxonomy(self, p1):
        _fetch_taxonomy(self.db_dir, threads=3, use_ftp=True)

        exp_cmd = [
            "kraken2-build", "--download-taxonomy",
            "--threads", "3", "--db", self.db_dir, "--use-ftp"
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
            _fetch_taxonomy(self.db_dir, threads=3, use_ftp=True)

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
            self.db_dir, libraries=libraries, all_kwargs=all_kwargs
        )

        base_cmd = ["kraken2-build", "--download-library"]
        exp_common_args = ["--threads", "2", "--db", self.db_dir]
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
                self.db_dir, libraries=['human'], all_kwargs=self.kwargs
            )

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_add_seqs_to_library(self, p1):
        seqs = DNAFASTAFormat(self.get_data_path('mags/samp1/bin1.fa'), 'r')

        _add_seqs_to_library(self.db_dir, seqs=seqs, no_masking=True)

        exp_cmd = [
            "kraken2-build", "--add-to-library", str(seqs.path),
            "--db", self.db_dir, "--no-mask"
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
            _add_seqs_to_library(self.db_dir, seqs=seqs, no_masking=True)

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_build_database(self, p1):
        _build_database(self.db_dir, all_kwargs=self.kwargs)

        exp_cmd = [
            "kraken2-build", "--build", "--db", self.db_dir,
            "--threads", "2", "--fast-build", "--kmer-len", "31",
            "--load-factor", "0.5", "--max-db-size", "1000"
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch('q2_moshpit.kraken2.database.run_command')
    def test_build_database_no_max_db(self, p1):
        all_kwargs = deepcopy(self.kwargs)
        all_kwargs['max_db_size'] = 0

        _build_database(self.db_dir, all_kwargs=all_kwargs)

        exp_cmd = [
            "kraken2-build", "--build", "--db", self.db_dir,
            "--threads", "2", "--fast-build", "--kmer-len", "31",
            "--load-factor", "0.5"
        ]
        p1.assert_called_once_with(cmd=exp_cmd, verbose=True)

    @patch(
        'q2_moshpit.kraken2.database.run_command',
        side_effect=CalledProcessError(123, 'cmd')
    )
    def test_build_database_exception(self, p1):
        with self.assertRaisesRegex(
                Exception,
                "An error was encountered .* building the database, "
                r"\(return code 123\), please inspect .*"
        ):
            _build_database(self.db_dir, all_kwargs=self.kwargs)

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

    @patch("q2_moshpit.kraken2.database.run_command")
    def test_build_kraken_db_action_standard_plus_library(self, p1):
        with self.assertRaisesRegex(
            ValueError,
            "Standard Kraken2 database was requested but some libraries"
        ):
            moshpit.actions.build_kraken_db(
                standard=True, libraries=['human'], threads=2
            )

    @patch("q2_moshpit.kraken2.database.Kraken2DBDirectoryFormat")
    @patch("q2_moshpit.kraken2.database._build_database")
    @patch("q2_moshpit.kraken2.database._add_seqs_to_library")
    @patch("q2_moshpit.kraken2.database._fetch_libraries")
    @patch("q2_moshpit.kraken2.database._fetch_taxonomy")
    @patch("q2_moshpit.kraken2.database._move_db_files")
    @patch("q2_moshpit.kraken2.database._build_standard_db")
    def test_build_kraken_db_action_with_existing_library(
            self, p1, p2, p3, p4, p5, p6, p7
    ):
        lib_path = 'path/to/lib'
        fake_dir_fmt = Kraken2DBDirectoryFormat(
            self.get_data_path('db'), 'r'
        )
        p7.return_value = fake_dir_fmt

        moshpit.actions.build_kraken_db(
            libraries=['human'], library_path=lib_path,
            threads=2, fast_build=True
        )

        p1.assert_not_called()
        p3.assert_called_once_with(
            db_dir=lib_path, threads=2, use_ftp=False
        )
        p4.assert_called_once_with(
            db_dir=lib_path, libraries=['human'], all_kwargs=ANY
        )
        p5.assert_not_called()
        p6.assert_called_once_with(
            db_dir=lib_path, all_kwargs=ANY
        )
        p2.assert_called_once_with(
            source=lib_path, destination=str(fake_dir_fmt.path)
        )

    @patch("q2_moshpit.kraken2.database.Kraken2DBDirectoryFormat")
    @patch("q2_moshpit.kraken2.database._build_database")
    @patch("q2_moshpit.kraken2.database._add_seqs_to_library")
    @patch("q2_moshpit.kraken2.database._fetch_libraries")
    @patch("q2_moshpit.kraken2.database._fetch_taxonomy")
    @patch("q2_moshpit.kraken2.database._move_db_files")
    @patch("q2_moshpit.kraken2.database._build_standard_db")
    def test_build_kraken_db_action_with_standard_library(
            self, p1, p2, p3, p4, p5, p6, p7
    ):
        lib_path = 'path/to/lib'
        fake_dir_fmt = Kraken2DBDirectoryFormat(
            self.get_data_path('db'), 'r'
        )
        p7.return_value = fake_dir_fmt

        moshpit.actions.build_kraken_db(
            standard=True, library_path=lib_path,
            threads=2, fast_build=True
        )

        p1.assert_called_once_with(
            db_dir=lib_path, all_kwargs=ANY
        )
        p2.assert_called_once_with(
            source=lib_path, destination=str(fake_dir_fmt.path)
        )
        p3.assert_not_called()
        p4.assert_not_called()
        p5.assert_not_called()
        p6.assert_not_called()

    @patch("q2_moshpit.kraken2.database.tempfile.TemporaryDirectory")
    @patch("q2_moshpit.kraken2.database.Kraken2DBDirectoryFormat")
    @patch("q2_moshpit.kraken2.database._build_database")
    @patch("q2_moshpit.kraken2.database._add_seqs_to_library")
    @patch("q2_moshpit.kraken2.database._fetch_libraries")
    @patch("q2_moshpit.kraken2.database._fetch_taxonomy")
    @patch("q2_moshpit.kraken2.database._move_db_files")
    @patch("q2_moshpit.kraken2.database._build_standard_db")
    def test_build_kraken_db_action_with_standard_library_temp(
            self, p1, p2, p3, p4, p5, p6, p7, p8
    ):
        fake_dir_fmt = Kraken2DBDirectoryFormat(
            self.get_data_path('db'), 'r'
        )
        p7.return_value = fake_dir_fmt
        fake_tmp_dir = MockTempDir()
        p8.return_value = fake_tmp_dir

        moshpit.actions.build_kraken_db(
            standard=True, threads=2, fast_build=True
        )

        p1.assert_called_once_with(
            db_dir=str(fake_tmp_dir.name), all_kwargs=ANY
        )
        p2.assert_called_once_with(
            source=str(fake_tmp_dir.name), destination=str(fake_dir_fmt.path)
        )
        p3.assert_not_called()
        p4.assert_not_called()
        p5.assert_not_called()
        p6.assert_not_called()

    @patch("q2_moshpit.kraken2.database.Kraken2DBDirectoryFormat")
    @patch("q2_moshpit.kraken2.database._build_database")
    @patch("q2_moshpit.kraken2.database._add_seqs_to_library")
    @patch("q2_moshpit.kraken2.database._fetch_libraries")
    @patch("q2_moshpit.kraken2.database._fetch_taxonomy")
    @patch("q2_moshpit.kraken2.database._move_db_files")
    @patch("q2_moshpit.kraken2.database._build_standard_db")
    def test_build_kraken_db_action_with_more_seqs(
            self, p1, p2, p3, p4, p5, p6, p7
    ):
        lib_path = 'path/to/lib'
        seqs = Artifact.import_data(
            'FeatureData[Sequence]', self.get_data_path("seqs")
        )
        fake_dir_fmt = Kraken2DBDirectoryFormat(
            self.get_data_path('db'), 'r'
        )
        p7.return_value = fake_dir_fmt

        moshpit.actions.build_kraken_db(
            seqs=[seqs], libraries=['human'], library_path=lib_path,
            threads=2, fast_build=True
        )

        p1.assert_not_called()
        p3.assert_called_once_with(
            db_dir=lib_path, threads=2, use_ftp=False
        )
        p4.assert_called_once_with(
            db_dir=lib_path, libraries=['human'], all_kwargs=ANY
        )
        p5.assert_called_once_with(
            db_dir=lib_path, seqs=ANY, no_masking=False
        )
        p6.assert_called_once_with(
            db_dir=lib_path, all_kwargs=ANY
        )
        p2.assert_called_once_with(
            source=lib_path, destination=str(fake_dir_fmt.path)
        )


if __name__ == "__main__":
    unittest.main()
