# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------



from q2_moshpit.eggnog import create_reference_db
from q2_types_genomics.feature_data import DiamondDB, MMseq2DB
from  q2_types.feature_data import FeatureData

from qiime2.plugin.testing import TestPluginBase

class TestEggnogDatabases(TestPluginBase):
    package = "q2_moshpit.tests"

    # TODO implement underlying functionality
    def test_detects_existing(self):
        pass


    def test_diamond_created(self):
        # test that output is of FeatureData[DiamondDB]
        pass


    def test_MMseq2_created(self):
        # test that output is of FeatureData[MMseq2DB]
        pass

class TestCreateDB(TestPluginBase):
    package = 'q2_moshpit.eggnog'

    def test_db_downloader_registered(self):
        #with self.assertRaises():
        create_reference_db(mode='diamond', target_taxa="bacteria",
                name="download_tester")

    def test_raises_taxa_parsing_error_with_taxids(self):
        # arrange
        bad_taxa_string = "72274,1123487" # TODO switch taxa values so they
        # will not work/raise the value error in `_check_taxa`

        # act
        # no further action needed, providing bad data and asserting as below
        # only things required.

        # assert
        with self.assertRaisesRegex(ValueError, "Unable to parse provided"
                " `Taxa` values"):
            create_reference_db("bad_taxa_test", )


