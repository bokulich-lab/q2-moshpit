# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import filecmp
from qiime2.plugin.testing import TestPluginBase
from q2_types.genome_data import OrthologAnnotationDirFmt
from ..annotations import collate_annotations


class TestCollateAnnotations(TestPluginBase):
    package = 'q2_moshpit.partition.tests'

    def test_collate_annotations(self):
        p = self.get_data_path("annotations")
        annotations = [
          OrthologAnnotationDirFmt(f"{p}/{letter}", mode="r")
          for letter in ["a", "b", "c"]
        ]
        collated_annotations = collate_annotations(annotations)

        # assert that all files are there
        compare = filecmp.dircmp(
            collated_annotations.path,
            self.get_data_path("annotations/collated")
        )
        self.assertListEqual(
            compare.common,
            [f"{letter}.annotations" for letter in ["a", "b", "c"]]
        )
