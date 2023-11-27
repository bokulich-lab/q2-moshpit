# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from ._method import (eggnog_diamond_search, eggnog_annotate, build_diamond_db)

__all__ = ['eggnog_diamond_search', 'eggnog_annotate', 'build_diamond_db']
