# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from collections import defaultdict
from urllib.parse import quote

import biom
import numpy as np
import pandas as pd
import scipy.sparse as sp

from .types import TFA_TAXONOMY_COLUMNS


def _gene_taxonomy_id(taxon_id: str, gene_id: str) -> str:
    """Join readable IDs, escaping the separator and literal percent signs."""
    return f"{quote(str(taxon_id), safe='')}|{quote(str(gene_id), safe='')}"


def _map_contigs_to_taxa(
    taxonomy: pd.DataFrame, taxon_to_contig_map: dict
) -> dict[str, str]:
    """Invert assignments for taxa with taxonomy, rejecting ambiguous contigs."""
    contig_to_taxon = {}
    for taxon_id, contigs in taxon_to_contig_map.items():
        if taxon_id not in taxonomy.index:
            continue
        for contig_id in contigs:
            if contig_id in contig_to_taxon and contig_to_taxon[contig_id] != taxon_id:
                raise ValueError(f"Contig {contig_id!r} is assigned to multiple taxa.")
            contig_to_taxon[contig_id] = taxon_id
    return contig_to_taxon


def _align_contig_matrices(
    contig_abundance: biom.Table,
    feature_inventory: biom.Table,
    contig_to_taxon: dict[str, str],
) -> tuple[list[str], sp.csr_matrix, sp.csr_matrix]:
    """Align shared contigs as rows of abundance and gene inventory matrices."""
    abundance_ids = list(contig_abundance.ids(axis="observation"))
    inventory_ids = list(feature_inventory.ids(axis="sample"))
    contig_ids = sorted(set(abundance_ids) & set(inventory_ids) & set(contig_to_taxon))
    abundance_index = {identifier: row for row, identifier in enumerate(abundance_ids)}
    inventory_index = {
        identifier: column for column, identifier in enumerate(inventory_ids)
    }

    # Both matrices have the same contig rows; columns are samples and genes.
    abundances = contig_abundance.matrix_data.tocsr()[
        [abundance_index[identifier] for identifier in contig_ids], :
    ]
    inventory = feature_inventory.matrix_data.tocsc()[
        :, [inventory_index[identifier] for identifier in contig_ids]
    ].T.tocsr()
    inventory.eliminate_zeros()
    return contig_ids, abundances, inventory


def _calculate_taxon_gene_loads(
    abundances: sp.csr_matrix,
    inventory: sp.csr_matrix,
    contig_ids: list[str],
    contig_to_taxon: dict[str, str],
    gene_ids: list[str],
) -> tuple[sp.csr_matrix, list[tuple[str, str]]]:
    """Calculate pair-by-sample loads, omitting pairs with no sample load."""
    contig_rows_by_taxon = defaultdict(list)
    for row, contig_id in enumerate(contig_ids):
        contig_rows_by_taxon[contig_to_taxon[contig_id]].append(row)

    blocks, pairs = [], []
    for taxon_id, contig_rows in sorted(contig_rows_by_taxon.items()):
        taxon_inventory = inventory[contig_rows, :]
        gene_columns = np.unique(taxon_inventory.indices)
        if not len(gene_columns):
            continue
        # Sparse multiplication avoids a dense taxon x gene x sample cube.
        loads = (
            taxon_inventory[:, gene_columns].T @ abundances[contig_rows, :]
        ).tocsr()
        loads.eliminate_zeros()
        nonzero_rows = np.flatnonzero(np.diff(loads.indptr))
        if len(nonzero_rows):
            blocks.append(loads[nonzero_rows, :])
            pairs.extend(
                (taxon_id, gene_ids[gene_columns[row]]) for row in nonzero_rows
            )

    matrix = (
        sp.vstack(blocks, format="csr")
        if blocks
        else sp.csr_matrix((0, abundances.shape[1]))
    )
    if not np.all(np.isfinite(matrix.data)):
        raise ValueError("Estimate loads must contain non-negative floats.")
    return matrix, pairs


def _estimate_tfa(
    contig_abundance: biom.Table,
    feature_inventory: biom.Table,
    taxonomy: pd.DataFrame,
    taxon_to_contig_map: dict,
) -> (biom.Table, pd.DataFrame):
    """Estimate sparse, sample-resolved taxon/function loads."""
    contig_to_taxon = _map_contigs_to_taxa(taxonomy, taxon_to_contig_map)
    contig_ids, abundances, inventory = _align_contig_matrices(
        contig_abundance, feature_inventory, contig_to_taxon
    )
    loads, pairs = _calculate_taxon_gene_loads(
        abundances,
        inventory,
        contig_ids,
        contig_to_taxon,
        list(feature_inventory.ids(axis="observation")),
    )

    if not pairs:
        raise ValueError("No taxon/gene pairs have nonzero loads in the input samples.")

    pair_ids = [_gene_taxonomy_id(taxon_id, gene_id) for taxon_id, gene_id in pairs]
    gene_taxonomy = pd.DataFrame(
        [
            [taxonomy.at[taxon_id, "Taxon"], taxon_id, gene_id]
            for taxon_id, gene_id in pairs
        ],
        index=pd.Index(pair_ids, name="Feature ID"),
        columns=TFA_TAXONOMY_COLUMNS,
    )
    return (
        biom.Table(
            loads,
            observation_ids=pair_ids,
            sample_ids=contig_abundance.ids(axis="sample"),
        ),
        gene_taxonomy,
    )


def estimate_tfa(
    ctx,
    contig_abundance,
    feature_inventory,
    taxonomy,
    taxon_to_contig_map,
    metadata=None,
):
    """Estimate sample-resolved TFA and visualize it with optional metadata."""
    estimate = ctx.get_action("annotate", "_estimate_tfa")
    explore = ctx.get_action("annotate", "_explore_tfa")
    tfa, gene_taxonomy = estimate(
        contig_abundance=contig_abundance,
        feature_inventory=feature_inventory,
        taxonomy=taxonomy,
        taxon_to_contig_map=taxon_to_contig_map,
    )
    (visualization,) = explore(
        feature_load=tfa, gene_taxonomy=gene_taxonomy, metadata=metadata
    )
    return tfa, gene_taxonomy, visualization
