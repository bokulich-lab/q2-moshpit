import numpy as np
import pandas as pd
from biom import Table
from rachis import Artifact


def factory_contig_abundance():
    return Artifact.import_data(
        "FeatureTable[Frequency]",
        Table(
            np.array(
                [
                    [100, 0, 0],
                    [50, 0, 0],
                    [0, 120, 0],
                    [0, 40, 0],
                    [0, 0, 90],
                    [0, 0, 60],
                ]
            ),
            ["contig1", "contig2", "contig3", "contig4", "contig5", "contig6"],
            ["sample1", "sample2", "sample3"],
        ),
    )


def factory_feature_inventory():
    return Artifact.import_data(
        "FeatureTable[Frequency]",
        Table(
            np.array([[2, 1, 2, 1, 2, 1], [1, 0, 1, 0, 1, 0], [0, 3, 0, 3, 0, 3]]),
            ["bla_TEM", "vanA", "mecA"],
            ["contig1", "contig2", "contig3", "contig4", "contig5", "contig6"],
        ),
    )


def factory_taxonomy():
    return Artifact.import_data(
        "FeatureData[Taxonomy]",
        pd.DataFrame(
            {
                "Taxon": [
                    "d__Bacteria;p__Actinobacteria;c__Actinomycetia;"
                    "o__Corynebacteriales;f__Mycobacteriaceae;"
                    "g__Mycobacterium;s__Mycobacterium tuberculosis",
                    "d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;"
                    "o__Enterobacteriales;f__Enterobacteriaceae;"
                    "g__Klebsiella;s__Klebsiella pneumoniae",
                ]
            },
            index=pd.Index(["T1", "T2"], name="Feature ID"),
        ),
    )


def factory_taxonomy_map():
    return Artifact.import_data(
        "FeatureMap[TaxonomyToContigs]",
        {
            "T1": ["contig1", "contig3", "contig5"],
            "T2": ["contig2", "contig4", "contig6"],
        },
    )


def estimate_tfa_example(use):
    contig_abundance = use.init_artifact("contig_abundance", factory_contig_abundance)
    feature_inventory = use.init_artifact(
        "feature_inventory", factory_feature_inventory
    )
    taxonomy = use.init_artifact("taxonomy", factory_taxonomy)
    taxon_to_contig_map = use.init_artifact("taxon_to_contig_map", factory_taxonomy_map)

    tfa, gene_taxonomy, visualization = use.action(
        use.UsageAction(plugin_id="annotate", action_id="estimate_tfa"),
        use.UsageInputs(
            contig_abundance=contig_abundance,
            feature_inventory=feature_inventory,
            taxonomy=taxonomy,
            taxon_to_contig_map=taxon_to_contig_map,
        ),
        use.UsageOutputNames(
            tfa="tfa", gene_taxonomy="gene_taxonomy", visualization="visualization"
        ),
    )

    tfa.assert_output_type("FeatureTable[Frequency % Properties('tfa')]")
    gene_taxonomy.assert_output_type("FeatureData[Taxonomy % Properties('tfa')]")
    visualization.assert_output_type("Visualization")
