# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
from collections import deque

from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat
)

import pandas as pd
import skbio


def kraken2_to_mag_features(kraken_reports: Kraken2ReportDirectoryFormat,
                           kraken_outputs: Kraken2OutputDirectoryFormat,
                           coverage_threshold: float = 0.1) \
         -> (pd.DataFrame, pd.DataFrame):
    table, taxonomy = kraken2_to_features(
        kraken_reports, coverage_threshold)

    rows_list = []
    taxa_list = []
    # convert IDs to match MAGs instead of taxids/db ids
    for sample_id in table.index:
        kraken_table_fp = (
            kraken_outputs.path / sample_id / f'{sample_id}.output.txt')
        hits_df = pd.read_csv(kraken_table_fp, sep='\t',
                              header=None, dtype='str')
        MAG_COL = 1
        TAXA_COL = 2

        sample_series = table.loc[sample_id, :]
        sample_obs = sample_series[sample_series != 0]
        merged_df = hits_df.join(sample_obs, on=TAXA_COL, how='right')
        merged_df = merged_df.join(taxonomy, on=TAXA_COL, how='left')

        new_taxa = merged_df[[MAG_COL, 'Taxon']].set_index(MAG_COL)
        new_taxa.index.name = 'Feature ID'

        table_row = pd.Series(True, index=new_taxa.index.unique())
        table_row.name = sample_id

        taxa_list.append(new_taxa)
        rows_list.append(table_row)

    cat_taxonomy = pd.concat(taxa_list)
    mag_taxonomy = cat_taxonomy[~cat_taxonomy.index.duplicated()]

    mag_table = pd.DataFrame(rows_list).fillna(False)

    return mag_table, mag_taxonomy


def kraken2_to_features(kraken_reports: Kraken2ReportDirectoryFormat,
                        coverage_threshold: float = 0.1) \
        -> (pd.DataFrame, pd.DataFrame):

    rows = []
    trees = []
    for relpath, df in kraken_reports.reports.iter_views(pd.DataFrame):
        sample_id, _ = os.path.split(relpath)

        filtered = df[df['perc_frags_covered'] >= coverage_threshold]
        tree = _kraken_to_ncbi_tree(filtered)
        tips = _ncbi_tree_to_tips(tree)
        if tips:
            table_row = pd.Series(True, index=_ncbi_tree_to_tips(tree))
            table_row.name = sample_id
            rows.append(table_row)
        trees.append(tree)

    full_tree = _combine_ncbi_trees(trees)

    table = pd.DataFrame(rows).fillna(False)
    taxonomy = _to_taxonomy(full_tree)

    return table, taxonomy


def _get_indentation(string, indent=2):
    return (len(string) - len(string.lstrip(' '))) // indent


def _kraken_to_ncbi_tree(df):
    tree = skbio.TreeNode()
    stack = deque([(0, tree)])
    for _, row in df.iterrows():
        r = row['rank']
        label = row['name']
        otu = str(row['ncbi_tax_id'])

        if r in ('U', 'R'):
            continue  # unclassified or root

        indent = _get_indentation(label)
        name = f"{r.lower()}__{label.strip()}"
        node = skbio.TreeNode(name=name, length=0)

        # Don't include internal non-strain infra-clades as tips
        if len(r) == 1 or r.startswith('S'):
            id_node = skbio.TreeNode(name=otu, length=0)
            node.length = 1  # not infra-clade, so give it a length
            node.append(id_node)

        parent_indent, parent_node = stack[-1]
        if parent_indent >= indent and parent_node.children:
            parent_node.children[0].is_actual_tip = True

        while parent_indent >= indent:
            stack.pop()
            parent_indent, parent_node = stack[-1]

        parent_node.append(node)
        stack.append((indent, node))

    # last entry is always a tip
    _, parent_node = stack[-1]
    if parent_node.children:
        parent_node.children[0].is_actual_tip = True

    return tree


def _combine_ncbi_trees(trees):
    full_tree = trees[0]
    for tree in trees[1:]:
        for tip in list(tree.tips()):
            try:
                # check if taxid is already in this tree
                full_tree.find(tip.name)
                continue  # for clarity
            except skbio.tree.MissingNodeError:
                parents = list(tip.ancestors())[:-1]  # ignore unnamed root
                matching = full_tree
                while parents:
                    node = parents.pop()
                    try:
                        matching = matching.find(node.name)
                    except skbio.tree.MissingNodeError:
                        matching.append(node)
                        break
    return full_tree


def _ncbi_tree_to_tips(tree):
    return [n.name for n in tree.tips() if hasattr(n, 'is_actual_tip')]


def _pad_ranks(ranks):
    order = ['d', 'k', 'p', 'c', 'o', 'f', 'g', 's']
    available = {}
    taxonomy = []

    for rank in reversed(ranks):
        r, label = rank.split('__', 1)
        available[r] = label
        if len(r) > 1 and r.startswith('s'):
            taxonomy.append(rank)

    for r in reversed(order):
        label = available.get(r)

        if label is not None:
            taxonomy.append(f'{r}__{label}')
            last_good_label = f'{r}__{label}'
        elif taxonomy:
            if r == 'k' and available.get('d') in ('Bacteria', 'Archaea'):
                # it is too strange to propogate upwards for these 'kingdoms'
                taxonomy.append(f"k__{available['d']}")
            else:
                # smear the most specific label we have upwards
                taxonomy.append(f'{r}__containing {last_good_label}')

    return ';'.join(reversed(taxonomy))


def _to_taxonomy(tree):
    rows = [(node.name, _pad_ranks(ranks))
            for node, ranks in tree.to_taxonomy()]
    taxonomy = pd.DataFrame(rows, columns=['Feature ID', 'Taxon'])
    taxonomy = taxonomy.set_index('Feature ID')

    return taxonomy
