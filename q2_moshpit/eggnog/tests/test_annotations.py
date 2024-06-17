# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import pandas as pd
from qiime2.plugin.testing import TestPluginBase

from q2_moshpit.eggnog import extract_annotations
from q2_types.feature_data_mag import OrthologAnnotationDirFmt


class TestAnnotationExtraction(TestPluginBase):
    package = 'q2_moshpit.eggnog.tests'

    def setUp(self):
        super().setUp()
        self.eggnog_ftf = pd.DataFrame(
            data={
                "ortholog1": [2, 0, 1, 10],
                "ortholog2": [0, 0, 1, 7],
                "ortholog3": [5, 2, 1, 0],
                "ortholog4": [1, 0, 1, 0],
                "ortholog5": [1, 3, 0, 9],
            },
            index=[
                "b9b4ab71-8e5f-48d7-bb23-df2726df1393",
                "62e07985-2556-435c-9e02-e7f94b8df07d",
                "1e9ffc02-0847-4f2c-b1e2-3965a4a78b15",
                "ab4f5ff0-45a1-41c9-9711-620765d5e92c"
            ]
        )
        self.mags_tpm = pd.DataFrame(
            data={
                "b9b4ab71-8e5f-48d7-bb23-df2726df1393": [0, 150.0, 10.0],
                "62e07985-2556-435c-9e02-e7f94b8df07d": [20.5, 15.0, 0.01],
                "1e9ffc02-0847-4f2c-b1e2-3965a4a78b15": [10.0, 1.5, 20.0],
                "ab4f5ff0-45a1-41c9-9711-620765d5e92c": [0.5, 0, 15.0],
            },
            index=["sample1", "sample2", "sample3"]
        )
        self.annotations = OrthologAnnotationDirFmt(
            self.get_data_path("annotations/"), mode='r'
        )

    def test_extract_annotations(self):
        obs_ft = extract_annotations(
            ortholog_annotations=self.annotations,
            annotation="cog"
        )
        exp_ft = pd.DataFrame(
            data={
                "L": [2.0, 0.0, 10.0, 3.0],
                "F": [1.0, 2.0, 0.0, 5.0],
                "A": [1.0, 0.0, 7.0, 0.0]
            },
            index=pd.Index([
                "1e9ffc02-0847-4f2c-b1e2-3965a4a78b15",
                "62e07985-2556-435c-9e02-e7f94b8df07d",
                "ab4f5ff0-45a1-41c9-9711-620765d5e92c",
                "b9b4ab71-8e5f-48d7-bb23-df2726df1393"
            ], name="id")
        )
        pd.testing.assert_frame_equal(obs_ft, exp_ft)
