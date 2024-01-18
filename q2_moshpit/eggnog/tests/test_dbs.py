# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import tempfile
from unittest.mock import patch, call
from qiime2.plugin.testing import TestPluginBase
from .._dbs import (
    fetch_eggnog_db, build_custom_diamond_db, fetch_eggnog_proteins,
    fetch_diamond_db, build_eggnog_diamond_db, fetch_ncbi_taxonomy,
    _write_version_tsv
)
from q2_types.feature_data import ProteinSequencesDirectoryFormat
from q2_types_genomics.reference_db import (
    NCBITaxonomyDirFmt, EggnogProteinSequencesDirFmt, NCBITaxonomyVersionFormat
)


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


class TestBuildDiamondDB(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    @patch("subprocess.run")
    def test_build_custom_diamond_db_simple(self, subp_run):
        # Instantiate input
        sequences = ProteinSequencesDirectoryFormat()

        # Call function. Patching will make sure nothing is
        # actually ran
        diamond_db = build_custom_diamond_db(sequences)

        # Paths to inputs and outputs
        path_in = os.path.join(str(sequences), "protein-sequences.fasta")
        path_out = os.path.join(str(diamond_db), "ref_db.dmnd")

        # Check that command was called in the expected way
        cmd = [
            "diamond", "makedb",
            "--verbose", "--log",
            "--in", f"{path_in}",
            "--db", f"{path_out}",
            "--threads", "1",
            '--file-buffer-size', '67108864'
        ]

        # Check that commands is ran as expected
        subp_run.assert_called_once_with(cmd, check=True)

    @patch("subprocess.run")
    def test_build_custom_diamond_db_with_taxonomy(self, subp_run):
        # Instantiate input
        sequences = ProteinSequencesDirectoryFormat()
        taxonomy_data = NCBITaxonomyDirFmt()

        # Call function. Patching will make sure nothing is
        # actually ran
        diamond_db = build_custom_diamond_db(sequences, taxonomy_data)

        # Paths to inputs and outputs
        path_in = os.path.join(str(sequences), "protein-sequences.fasta")
        path_tax_map = os.path.join(
            str(taxonomy_data), "prot.accession2taxid.gz"
            )
        path_tax_nodes = os.path.join(str(taxonomy_data), "nodes.dmp")
        path_tax_names = os.path.join(str(taxonomy_data), "names.dmp")
        path_out = os.path.join(str(diamond_db), "ref_db.dmnd")

        # Check that command was called in the expected way
        cmd = [
            "diamond", "makedb",
            "--verbose", "--log",
            "--in", f"{path_in}",
            "--db", f"{path_out}",
            "--threads", "1",
            '--file-buffer-size', '67108864',
            "--taxonmap", f"{path_tax_map}",
            "--taxonnodes", f"{path_tax_nodes}",
            "--taxonnames", f"{path_tax_names}",
        ]

        # Check that commands is ran as expected
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

    @patch("q2_moshpit.eggnog._dbs._write_version_tsv")
    @patch("subprocess.run")
    def test_fetch_ncbi_taxonomy(self, subp_run, w_v_tsv):
        # Call function. Patching will make sure nothing is actually ran
        ncbi_data = fetch_ncbi_taxonomy()
        zip_path = os.path.join(str(ncbi_data), "taxdmp.zip")
        nodes_path = os.path.join(str(ncbi_data), "nodes.dmp")
        names_path = os.path.join(str(ncbi_data), "names.dmp")
        proteins_path = os.path.join(str(ncbi_data), "prot.accession2taxid.gz")
        version_path = os.path.join(str(ncbi_data), "version.tsv")

        # Check that command was called in the expected way
        first_call = call(
            [
                "wget", "-O", zip_path,
                "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
            ],
            check=True
        )
        second_call = call(
            [
                "unzip", "-j", zip_path, "names.dmp", "nodes.dmp",
                "-d", str(ncbi_data)
            ],
            check=True,
        )
        third_call = call(
            ["rm", zip_path],
            check=True,
        )
        forth_call = call(
            [
                "wget", "-O", proteins_path,
                "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/"
                "prot.accession2taxid.gz"
            ],
            check=True,
        )

        # Check that commands are ran as expected
        subp_run.assert_has_calls(
            [first_call, second_call, third_call, forth_call],
            any_order=False
        )
        w_v_tsv.assert_called_once_with(
            nodes_path,
            names_path,
            proteins_path,
            version_path
        )

    def test_make_version_df(self):
        nodes = self.get_data_path('ncbi/nodes.dmp')
        names = self.get_data_path('ncbi/names.dmp')
        proteins = self.get_data_path('ncbi/prot.accession2taxid.gz')

        with tempfile.TemporaryDirectory() as tmp:
            version = os.path.join(tmp, 'version.tsv')
            _write_version_tsv(nodes, names, proteins, version)
            format = NCBITaxonomyVersionFormat(version, mode="r")
            format.validate()

    @patch("q2_moshpit.eggnog._dbs._validate_taxon_id")
    @patch("subprocess.run")
    @patch("shutil.move")
    def test_build_eggnog_diamond_db(self, shut_mv, subp_run, _val):
        # Instantiate input
        proteins_and_taxa = EggnogProteinSequencesDirFmt()

        # Call function. Patching will make sure nothing is
        # actually ran
        diamond_db = build_eggnog_diamond_db(proteins_and_taxa, taxon=2)

        # Check that command was called in the expected way
        cmd = [
            "create_dbs.py",
            "--data_dir", str(proteins_and_taxa),
            "--taxids", "2",
            "--dbname", "ref_db"
        ]

        # Check that subprocess.run is run as expected
        subp_run.assert_called_once_with(cmd, check=True)

        # Check that shutil.move is run as expected
        source_path = os.path.join(str(proteins_and_taxa), "ref_db.dmnd")
        destination_path = os.path.join(str(diamond_db), "ref_db.dmnd")
        shut_mv.assert_called_once_with(source_path, destination_path)

    def test_build_eggnog_diamond_db_invalid_taxon_id(self):
        # Init input data
        path_to_data = self.get_data_path('build_eggnog_diamond_db/')
        eggnog_proteins = EggnogProteinSequencesDirFmt(path_to_data, 'r')

        # Call function exception error since taxon 0 is invalid
        with self.assertRaisesRegex(
            ValueError,
            "'0' is not valid taxon ID. "
        ):
            _ = build_eggnog_diamond_db(eggnog_proteins, taxon=0)
