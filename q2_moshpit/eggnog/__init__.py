# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from ._method import (
    eggnog_diamond_search, _eggnog_diamond_search, eggnog_annotate,
    _eggnog_feature_table, _eggnog_annotate
)
from ._dbs import (
    fetch_eggnog_db, fetch_diamond_db, build_custom_diamond_db,
    fetch_eggnog_proteins, build_eggnog_diamond_db, fetch_ncbi_taxonomy
)


__all__ = [
    'eggnog_diamond_search', '_eggnog_diamond_search', 'eggnog_annotate',
    '_eggnog_feature_table', 'fetch_eggnog_db', 'fetch_diamond_db',
    'build_custom_diamond_db', 'fetch_eggnog_proteins',
    'build_eggnog_diamond_db', 'fetch_ncbi_taxonomy', '_eggnog_annotate'
]
