# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from .dereplication import dereplicate_mags
from .kaiju import classification as kaiju_class, database as kaiju_db
from .kraken2 import (
    classification as kraken_class, database as kraken_db, bracken, helpers
)
from .metabat2 import metabat2
from . import prodigal
from . import eggnog
from . import busco


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

__all__ = [
    'metabat2', 'bracken', 'kraken_class', 'kraken_db',
    'kaiju_class', 'kaiju_db', 'dereplicate_mags', 'eggnog',
    'busco', 'prodigal', 'helpers'
]
