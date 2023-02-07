# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


import importlib
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

importlib.import_module('q2_moshpit.kraken2')
importlib.import_module('q2_moshpit.metabat2')
importlib.import_module('q2_moshpit.diamond')
importlib.import_module('q2_moshpit.annotation')
