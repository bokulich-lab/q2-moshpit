# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from q2_moshpit.eggnog import e_mapper
from q2_types_genomics.feature_data import EggnogRefDirFmt
from q2_types.feature_data import DNAFASTAFormat
from qiime2.plugin.testing import TestPluginBase
from qiime2 import Artifact


class TestMapper(TestPluginBase):
    package = "q2_moshpit.eggnog.tests"

    def test_basic_passing(self):
        query_artifact = Artifact.load(
            self.get_data_path("query_sequences.qza"))
        query_seqs = query_artifact.view(DNAFASTAFormat)

        main_ref_artifact = Artifact.load(self.get_data_path(
            "refrence_dbs.qza"))
        main_ref = main_ref_artifact.view(EggnogRefDirFmt)

        e_mapper(query_sequences=query_seqs,
                 main_db=main_ref, itype='genome',
                 )
