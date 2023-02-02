# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .kraken2 import classification, database
from .metabat2 import metabat2

import importlib
from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

importlib.import_module('q2_moshpit.metabat2')
importlib.import_module('q2_moshpit.diamond')
importlib.import_module('q2_moshpit.annotation')
# importlib.import_module('q2_moshpit.usage_examples')
__all__ = ['metabat2', 'classification', 'database']
