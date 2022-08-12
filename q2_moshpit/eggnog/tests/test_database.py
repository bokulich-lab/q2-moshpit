# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------



from q2_moshpit.eggnog import create_reference_db, _check_taxa
from q2_types_genomics.feature_data import DiamondDB, MMseq2DB
from  q2_types.feature_data import FeatureData

from qiime2.plugin.testing import TestPluginBase

# set temp data directory up...


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

    def test_db_downloader_runs(self):
        exp = """Parsing taxa table ...
Total taxa selected: 0
Parsing fasta file ...
Total sequences selected: 0
Total sequences discarded: 0"""

        #with self.assertRaises():
        create_reference_db(mode='diamond', target_taxa="bacteria",
                name="download_tester", simulate=True)

    def test_raises_taxa_parsing_error_with_taxids(self):
        # arrange
        bad_taxa_string = "72274,1123487,wombat" # TODO switch taxa values so they
        # will not work/raise the value error in `_check_taxa`

        # act
        # no further action needed, providing bad data and asserting as below
        # only things required.

        # assert
        with self.assertRaisesRegex(ValueError, "Unable to parse provided"
                " `Taxa` values"):
            create_reference_db(mode='diamond', target_taxa=bad_taxa_string,
                    simulate=True)

class TestDButils(TestPluginBase):
    package = 'q2_moshpit.eggnog'
    def testraise_on_mixed_numeric(self):
        pass
    def test_raise_on_mixed_string(self):
        pass
    def test_successful_numeric(self):
        pass
    def test_successful_string(self):
        pass
