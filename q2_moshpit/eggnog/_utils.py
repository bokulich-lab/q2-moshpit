# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------


def _check_taxa(target_taxa: list):
    if target_taxa[0].isnumeric():
        taxa_type = "--taxids"
        for taxon in target_taxa:
            if not taxon.isnumeric():
                raise ValueError("All taxa inputs must be the same type,"
                                 " either all taxids as integers or all"
                                 " string labels.")

    else:
        taxa_type = "--taxa"
        for taxon in target_taxa:
            if taxon.isnumeric():
                raise ValueError(
                    "All taxa inputs must be the same type,"
                    " either all taxids as integers or all string labels.")

    return taxa_type, ",".join(target_taxa)
