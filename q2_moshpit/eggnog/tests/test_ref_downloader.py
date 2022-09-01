# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from q2_moshpit.eggnog import create_reference_db, download_references

from pytest import mark
from qiime2.plugin.testing import TestPluginBase

# set temp data directory up...


class TestEggnogDatabases(TestPluginBase):
    package = "q2_moshpit.eggnog.tests"

    # TODO implement underlying functionality
    def test_detects_existing(self):
        pass

    def test_diamond_created(self):
        # test that output is of FeatureData[DiamondDB]
        pass

    def test_MMseq2_created(self):
        # test that output is of FeatureData[MMseq2DB]
        pass


class TestRefDownload(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    def test_all_download(self):
        download_references()


class TestCreateDB(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    def test_db_downloader_diamond_basic(self):
        # with self.assertRaises():
        create_reference_db(mode='diamond', target_taxa="bacteria",
                            name="download_tester", simulate=True)

    def test_raises_taxa_parsing_error_with_taxids(self):
        bad_taxa_string = ["72274", "1123487", "wombat"]

        with self.assertRaisesRegex(ValueError, "All taxa inputs must be the"
                                    " same type, either all taxids as"
                                    " integers or all string labels."):
            create_reference_db(mode='diamond', target_taxa=bad_taxa_string,
                                simulate=True)

    def test_diamond_passes_all_numeric(self):
        numeric_taxa = ['123', '879', '5598777']
        create_reference_db(mode='diamond', target_taxa=numeric_taxa,
                            simulate=True)

    def test_diamond_passes_all_named(self):
        named_taxa = ['Gammaproteobacteria', 'Aquificae', 'Persephonella']
        create_reference_db(mode='diamond', target_taxa=named_taxa,
                            simulate=True)

    def test_bad_mode_failing(self):
        with self.assertRaisesRegex(ValueError, "Please supply a valid mode"):
            create_reference_db(mode='pancake', target_taxa="bacteria",
                                name="bad_mode", simulate=True)

    @mark.skip(reason="stupendously slow for some reason....")
    def test_actual_download_diamond(self):
        create_reference_db(mode="diamond", target_taxa="85004",
                            name="real_download")

    @mark.skip(reason="stupendously slow for some reason....")
    def test_actual_download_mmseqs(self):
        create_reference_db(mode="mmseqs", target_taxa="85004",
                            name="real_download")
