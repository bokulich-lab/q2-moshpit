# ----------------------------------------------------------------------------
# Copyright (c) 2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.plugin.testing import TestPluginBase
from q2_moshpit.eggnog._format import EggnogHmmerIdmapFileFmt
from qiime2.plugin import ValidationError


class TestEggnogHmmerIdmap(TestPluginBase):
    package = "q2_moshpit.eggnog.tests"

    def test_HmmIdmapFileFmt_valid(self):
        fmt = EggnogHmmerIdmapFileFmt(
            self.get_data_path("valid_idmaps/bacteria.hmm.idmap"), 'r'
        )
        fmt.validate()

    def test_HmmIdmapFileFmt_invalid_idmap_1(self):
        fmt = EggnogHmmerIdmapFileFmt(
            self.get_data_path("invalid_idmaps/1.hmm.idmap"), 'r'
        )
        with self.assertRaisesRegex(
            ValidationError,
            "Expected index and an alphanumeric code separated "
            "by a single space."
        ):
            fmt.validate(level="min")

    def test_HmmIdmapFileFmt_invalid_idmap_2(self):
        fmt = EggnogHmmerIdmapFileFmt(
            self.get_data_path("invalid_idmaps/2.hmm.idmap"), 'r'
        )
        with self.assertRaisesRegex(
            ValidationError,
            "Expected index and an alphanumeric code separated "
            "by a single space."
        ):
            fmt.validate(level="min")

    def test_HmmIdmapFileFmt_invalid_idmap_3(self):
        fmt = EggnogHmmerIdmapFileFmt(
            self.get_data_path("invalid_idmaps/3.hmm.idmap"), 'r'
        )
        with self.assertRaisesRegex(
            ValidationError,
            'Expected index'
        ):
            fmt.validate(level="min")

    def test_HmmIdmapFileFmt_invalid_idmap_4(self):
        fmt = EggnogHmmerIdmapFileFmt(
            self.get_data_path("invalid_idmaps/4.hmm.idmap"), 'r'
        )
        with self.assertRaisesRegex(
            ValidationError,
            "Expected index and an alphanumeric code separated "
            "by a single space."
        ):
            fmt.validate(level="min")
