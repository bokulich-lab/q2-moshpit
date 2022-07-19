# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


from ..plugin_setup import plugin
from qiime2.plugin import SemanticType
from q2_types_genomics.eggnog import ArbitraryHeaderTSVDirFmt


DiamondDB = SemanticType('DiamondDB')

MMseq2DB = SemanticType('MMseq2DB')

plugin.register_semantic_types(
        DiamondDB, MMseq2DB,
)

plugin.register_semantic_type_to_format(
    DiamondDB,
    artifact_format=ArbitraryHeaderTSVDirFmt
)

plugin.register_semantic_type_to_format(
    MMseq2DB,
    artifact_format=ArbitraryHeaderTSVDirFmt
)
