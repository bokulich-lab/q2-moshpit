# ----------------------------------------------------------------------------
# Copyright (c) 2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import json
from pathlib import Path
import qiime2
import biom
import numpy as np
import pandas as pd
import scipy.sparse as sp
from qiime2.plugin.testing import TestPluginBase

from ..tfa import (
    _align_contig_matrices,
    _gene_taxonomy_id,
    _map_contigs_to_taxa,
    _estimate_tfa,
)


class TestTFA(TestPluginBase):
    package = "q2_annotate.tfa.tests"

    def setUp(self):
        super().setUp()
        self.abundance = pd.read_csv(
            self.get_data_path("contig-abundance.tsv"), sep="\t", index_col=0
        ).T
        self.inventory = pd.read_csv(
            self.get_data_path("gene-inventory.tsv"), sep="\t", index_col=0
        ).T
        self.taxonomy = pd.DataFrame(
            {"Taxon": ["Escherichia coli", "Klebsiella pneumoniae"]},
            index=pd.Index(["T1", "T2"], name="Feature ID"),
        )
        self.mapping = json.load(open(self.get_data_path("taxon-to-contigs.json")))
        self.expected_tfa = pd.read_csv(
            self.get_data_path("expected-tfa-table.tsv"), sep="\t", index_col=0
        ).astype(float)
        self.expected_gene_taxonomy = pd.read_csv(
            self.get_data_path("expected-gene-taxonomy.tsv"),
            sep="\t",
            index_col=0,
            dtype=str,
            keep_default_na=False,
        )

    def _estimate(self, abundance=None, inventory=None, taxonomy=None, mapping=None):
        abundance = self.abundance if abundance is None else abundance
        inventory = self.inventory if inventory is None else inventory
        taxonomy = self.taxonomy if taxonomy is None else taxonomy
        mapping = self.mapping if mapping is None else mapping
        abundance, inventory = abundance.T, inventory.T
        return _estimate_tfa(
            biom.Table(
                abundance.values,
                observation_ids=abundance.index,
                sample_ids=abundance.columns,
            ),
            biom.Table(
                inventory.values,
                observation_ids=inventory.index,
                sample_ids=inventory.columns,
            ),
            taxonomy,
            mapping,
        )

    def test_sample_resolved_values_and_sparse_pairs(self):
        """Estimate per-sample loads while storing only observed taxon-gene pairs."""
        observed, gene_taxonomy = self._estimate()
        self.assertTrue(sp.issparse(observed.matrix_data))
        self.assertEqual(
            observed.matrix_data.nnz, np.count_nonzero(self.expected_tfa.values)
        )
        pd.testing.assert_frame_equal(
            observed.to_dataframe(dense=True).rename_axis("Feature ID"),
            self.expected_tfa,
        )
        pd.testing.assert_frame_equal(gene_taxonomy, self.expected_gene_taxonomy)
        self.assertEqual(
            list(observed.ids(axis="observation")), list(gene_taxonomy.index)
        )

    def test_no_common_contigs_rejected(self):
        """Report when no mapped contigs overlap the input tables."""
        with self.assertRaisesRegex(
            ValueError, "No taxon/gene pairs have nonzero loads"
        ):
            self._estimate(mapping={"T1": ["not-a-contig"]})

    def test_ignores_extra_contigs(self):
        """Ignore mapped contigs absent from the gene inventory."""
        abundance = self.abundance.assign(C7=[1000, 0, 0])
        taxonomy = pd.DataFrame(
            {"Taxon": ["Escherichia coli", "Klebsiella pneumoniae", "Other"]},
            index=pd.Index(["T1", "T2", "T3"]),
        )
        mapping = {**self.mapping, "T3": ["C7", "C8"]}
        observed, gene_taxonomy = self._estimate(
            abundance=abundance, taxonomy=taxonomy, mapping=mapping
        )
        pd.testing.assert_frame_equal(
            observed.to_dataframe(dense=True).rename_axis("Feature ID"),
            self.expected_tfa,
        )
        pd.testing.assert_frame_equal(gene_taxonomy, self.expected_gene_taxonomy)

    def test_zero_abundance_rejected(self):
        """Report when all pairs have zero loads and no taxonomy can be returned."""
        with self.assertRaisesRegex(
            ValueError, "No taxon/gene pairs have nonzero loads"
        ):
            self._estimate(abundance=self.abundance * 0)

    def test_omits_pairs_with_zero_load_across_samples(self):
        """Keep loaded pairs while dropping an entirely absent taxon."""
        abundance = self.abundance.copy()
        abundance[["C2", "C4", "C6"]] = 0
        observed, gene_taxonomy = self._estimate(abundance=abundance)
        self.assertEqual(
            list(observed.ids(axis="observation")), ["T1|bla_TEM", "T1|vanA"]
        )
        self.assertEqual(list(gene_taxonomy.index), ["T1|bla_TEM", "T1|vanA"])

    def test_extreme_abundances(self):
        """Preserve very small and very large abundance values."""
        abundance = self.abundance * 0.0
        abundance.loc["S1", "C1"] = 1e-9
        abundance.loc["S1", "C2"] = 1e12
        observed, gene_taxonomy = self._estimate(abundance=abundance)
        loads = {
            tuple(gene_taxonomy.loc[pair_id, ["Taxon ID", "Gene ID"]]): values
            for pair_id, values in zip(
                observed.ids(axis="observation"), observed.matrix_data.toarray()
            )
        }
        np.testing.assert_allclose(
            loads[("T1", "bla_TEM")],
            [2e-9, 0, 0],
            rtol=1e-12,
            atol=0,
        )
        np.testing.assert_allclose(
            loads[("T2", "mecA")],
            [3e12, 0, 0],
            rtol=1e-12,
            atol=0,
        )

    def test_mapping_preserves_ambiguous_characters(self):
        """Preserve punctuation in taxon and gene IDs through estimation."""
        inventory = self.inventory.rename(columns={"bla_TEM": 'a,"b'})
        taxonomy = self.taxonomy.rename(index={"T1": 't,"1'})
        mapping = {'t,"1': self.mapping["T1"], "T2": self.mapping["T2"]}
        observed, gene_taxonomy = self._estimate(
            inventory=inventory, taxonomy=taxonomy, mapping=mapping
        )
        pairs = set(
            gene_taxonomy[["Taxon ID", "Gene ID"]].itertuples(index=False, name=None)
        )
        self.assertIn(('t,"1', 'a,"b'), pairs)

    def test_ambiguous_taxon_assignment_rejected(self):
        """Reject a contig assigned to more than one taxon."""
        mapping = {**self.mapping, "T2": ["C1", "C2", "C4", "C6"]}
        with self.assertRaisesRegex(ValueError, "multiple taxa"):
            self._estimate(mapping=mapping)

    def test_contig_mapping_skips_absent_taxonomy(self):
        """Skip mapped taxa missing from the taxonomy table."""
        self.assertEqual(
            _map_contigs_to_taxa(self.taxonomy.loc[["T1"]], self.mapping),
            {"C1": "T1", "C3": "T1", "C5": "T1"},
        )

    def test_align_contig_matrices_preserves_axes_and_values(self):
        """Align contigs across matrices without changing sample or gene values."""
        abundance = self.abundance.T.iloc[::-1]
        inventory = self.inventory.T.iloc[:, ::-1]
        contigs, abundances, genes = _align_contig_matrices(
            biom.Table(
                abundance.values,
                observation_ids=abundance.index,
                sample_ids=abundance.columns,
            ),
            biom.Table(
                inventory.values,
                observation_ids=inventory.index,
                sample_ids=inventory.columns,
            ),
            {"C1": "T1", "C3": "T1", "C5": "T1"},
        )
        self.assertEqual(contigs, ["C1", "C3", "C5"])
        self.assertIsInstance(abundances, sp.csr_matrix)
        self.assertIsInstance(genes, sp.csr_matrix)
        np.testing.assert_array_equal(
            abundances.toarray(), self.abundance[contigs].T.values
        )
        np.testing.assert_array_equal(
            genes.toarray(), self.inventory.loc[contigs].values
        )

    def test_stable_feature_ids_and_mapping(self):
        """Keep feature IDs and loads stable when input axes are reordered."""
        table, mapping = self._estimate()
        reordered, reordered_mapping = self._estimate(
            inventory=self.inventory[self.inventory.columns[::-1]],
            abundance=self.abundance[self.abundance.columns[::-1]],
        )
        pd.testing.assert_frame_equal(
            mapping.sort_index(), reordered_mapping.sort_index()
        )
        self.assertEqual(list(table.ids(axis="observation")), list(mapping.index))
        self.assertEqual(set(mapping.index), set(self.expected_gene_taxonomy.index))
        pd.testing.assert_frame_equal(
            table.to_dataframe(dense=True).sort_index(),
            reordered.to_dataframe(dense=True).sort_index(),
        )

    def test_readable_ids_escape_ambiguous_components(self):
        """Build unique, readable feature IDs from taxon and gene IDs."""
        cases = pd.read_csv(self.get_data_path("pair-ids.tsv"), sep="\t", dtype=str)
        observed = [
            _gene_taxonomy_id(taxon, gene)
            for taxon, gene in zip(cases["Taxon ID"], cases["Gene ID"])
        ]
        self.assertEqual(observed, list(cases["Feature ID"]))
        self.assertEqual(len(observed), len(set(observed)))

    def _pipeline_inputs(self):
        """Import the same readable fixtures used for estimator tests."""
        abundance, inventory = self.abundance.T, self.inventory.T
        return {
            "contig_abundance": qiime2.Artifact.import_data(
                "FeatureTable[Frequency]",
                biom.Table(abundance.values, abundance.index, abundance.columns),
            ),
            "feature_inventory": qiime2.Artifact.import_data(
                "FeatureTable[Frequency]",
                biom.Table(inventory.values, inventory.index, inventory.columns),
            ),
            "taxonomy": qiime2.Artifact.import_data(
                "FeatureData[Taxonomy]", self.taxonomy
            ),
            "taxon_to_contig_map": qiime2.Artifact.import_data(
                "FeatureMap[TaxonomyToContigs]", self.mapping
            ),
        }

    def test_pipeline_registration(self):
        """Expose one pipeline and register its two private helper actions."""
        self.assertIn("estimate_tfa", self.plugin.pipelines)
        self.assertNotIn("estimate_tfa", self.plugin.methods)
        self.assertIn("_estimate_tfa", self.plugin.methods)
        self.assertIn("_explore_tfa", self.plugin.visualizers)
        self.assertNotIn("explore_tfa", self.plugin.visualizers)

    def test_pipeline_outputs_with_optional_metadata(self):
        """Return matching tables and an explorer with optional sample groups."""
        inputs = self._pipeline_inputs()
        metadata = qiime2.Metadata.load(self.get_data_path("sample-metadata.tsv"))
        for grouping in (None, metadata):
            with self.subTest(metadata=grouping is not None):
                result = self.plugin.pipelines["estimate_tfa"](
                    **inputs, metadata=grouping
                )
                self.assertEqual(len(result), 3)
                self.assertEqual(
                    str(result.tfa.type), "FeatureTable[Frequency % Properties('tfa')]"
                )
                self.assertEqual(
                    str(result.gene_taxonomy.type),
                    "FeatureData[Taxonomy % Properties('tfa')]",
                )
                pd.testing.assert_frame_equal(
                    result.tfa.view(biom.Table)
                    .to_dataframe(dense=True)
                    .rename_axis("Feature ID"),
                    self.expected_tfa,
                )
                pd.testing.assert_frame_equal(
                    result.gene_taxonomy.view(pd.DataFrame), self.expected_gene_taxonomy
                )
                self.assertIsInstance(result.visualization, qiime2.Visualization)
                result.visualization.export_data(self.temp_dir.name)
                html = (Path(self.temp_dir.name) / "index.html").read_text()
                payload = json.loads(
                    html.split('<script id="tfa-data" type="application/json">')[
                        1
                    ].split("</script>")[0]
                )
                self.assertEqual(payload["samples"], ["S1", "S2", "S3"])
                self.assertEqual(bool(payload["groups"]), grouping is not None)

    def test_pipeline_rejects_unobserved_pairs(self):
        """Propagate the estimator error when all sample loads are zero."""
        inputs = self._pipeline_inputs()
        table = inputs["contig_abundance"].view(biom.Table)
        inputs["contig_abundance"] = qiime2.Artifact.import_data(
            "FeatureTable[Frequency]",
            biom.Table(
                sp.csr_matrix(table.shape),
                table.ids(axis="observation"),
                table.ids(axis="sample"),
            ),
        )
        with self.assertRaisesRegex(
            ValueError, "No taxon/gene pairs have nonzero loads"
        ):
            self.plugin.pipelines["estimate_tfa"](**inputs)
