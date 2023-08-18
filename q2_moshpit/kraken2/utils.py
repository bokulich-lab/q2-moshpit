# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from collections import Counter
from itertools import zip_longest, takewhile
from re import sub

from q2_moshpit._utils import _construct_param


def _process_kraken2_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by Kraken 2.

    Argument names will be converted to command line parameters by
    appending a '--' prefix and replacing all '_' with '-',
    e.g.: 'some_parameter_x' -> '--some-parameter-x'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """
    if isinstance(arg_val, bool) and arg_val:
        return [_construct_param(arg_key)]
    elif not isinstance(arg_val, list):
        return [_construct_param(arg_key), str(arg_val)]
    else:
        raise NotImplementedError(
            f'Parsing arguments of type "{type(arg_val)}" is not supported.'
        )


def _find_lca(taxa):
    """Find least common ancestor between two semicolon-delimited strings."""
    # determine optimal zip mode. Normally, zip is best because trimming to
    # shortest is an inherent feature of LCA.
    # However if only one frame contains an assignment for feature x, we want
    # to just take that taxonomy. zip_longest will accomplish this while using
    # the same machinery... same deal if we want to take a superset taxonomy
    if '' in taxa:
        zip_it = zip_longest
    else:
        zip_it = zip
    # LCA ends where zipped taxonomy strings no longer converge to len == 1
    taxa_comparison = [set(rank) - {None} for rank in zip_it(*taxa)]
    return (rank.pop() for rank in takewhile(
        lambda x: len(x) == 1, taxa_comparison))


# LCA majority is same as super majority without substring collapsing
def _find_lca_majority(taxa):
    return _find_super_lca(taxa, collapse_substrings=False)


# modified version of _find_lca that prioritizes majority and supersets
def _find_super_lca(taxa, collapse_substrings=True):
    # collapse and count unique labels at each rank
    # yields list of ('labels', counts) sorted by most to least abundant
    taxa = [[t for t in r if t not in ['',]] for r in zip_longest(*taxa)]
    # taxa = [[t if t else '' for t in r] for r in taxa]
    if collapse_substrings:
        # find longest string in group of sub/superstrings, combine
        taxa = [[max({i for i in x if t in i}, key=len) for t in x]
                if x else '' for x in taxa]
    taxa_comparison = [Counter(t).most_common() for t in taxa]
    # return majority wherever a clear majority is found
    # terminate when no majority is found, that's your LCA
    # propagate empty ranks that maintain majority/consensus by inserting ''
    # to preserve empty ranks in the original taxonomies (e.g., when using a
    # rank handle unannotated levels will be removed, but as long as these pass
    # muster we want to preserve those labels)
    return [rank[0][0] if rank else '' for rank in takewhile(
        lambda x: len(x) < 2 or x[0][1] > x[1][1], taxa_comparison)]


def _taxon_to_list(taxon, rank_handle):
    """Split taxonomy string into list of taxonomic labels"""
    if rank_handle != '':
        return [sub(rank_handle, '', t.strip()) for t in taxon.split(';')]
    else:
        return [t.strip() for t in taxon.split(';')]
