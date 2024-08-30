# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from . import busco
from . import eggnog
from . import partition
from . import prodigal
from ._version import get_versions
from .dereplication import dereplicate_mags
from .filtering import filter_derep_mags, filter_mags, filter_reads_pangenome
from .kaiju import classification as kaiju_class, database as kaiju_db
from .kraken2 import (
    classification as kraken_class,
    database as kraken_db, bracken,
    helpers as kraken_helpers
)
from .metabat2 import metabat2
from ._utils import get_feature_lengths

__version__ = get_versions()['version']
del get_versions

__all__ = [
    'metabat2', 'bracken', 'kraken_class', 'kraken_db',
    'kaiju_class', 'kaiju_db', 'dereplicate_mags', 'eggnog',
    'busco', 'prodigal', 'kraken_helpers', 'partition',
    'filter_derep_mags', 'filter_mags', 'get_feature_lengths',
    'filter_reads_pangenome'
]
