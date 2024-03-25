# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import csv

from qiime2.core.exceptions import ValidationError
from qiime2.plugin import model


class BUSCOResultsFormat(model.TextFileFormat):
    HEADER = [
        "mag_id", "sample_id", "input_file", "dataset", "complete",
        "single", "duplicated", "fragmented", "missing", "n_markers",
        "scaffold_n50", "contigs_n50", "percent_gaps", "scaffolds",
        "length"
    ]

    def _validate(self, n_records=None):
        with self.open() as fh:
            reader = csv.reader(fh, delimiter='\t')
            headers = next(reader)

            if set(headers) != set(self.HEADER):
                raise ValidationError(
                    f'Invalid header: {headers}, expected: {self.HEADER}'
                )

            for i, row in enumerate(reader, start=2):
                if len(row) != len(self.HEADER):
                    raise ValidationError(f'Line {i} has {len(row)} columns, expected {len(self.HEADER)}')

                if n_records is not None and i - 1 >= n_records:
                    break

    def _validate_(self, level):
        record_count_map = {'min': 100, 'max': None}
        self._validate(record_count_map[level])


BUSCOResultsDirectoryFormat = model.SingleFileDirectoryFormat(
    'BUSCOResultsDirectoryFormat', 'busco_results.tsv',
    BUSCOResultsFormat
)
