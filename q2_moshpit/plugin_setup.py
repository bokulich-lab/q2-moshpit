# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

from q2_types.distance_matrix import DistanceMatrix
from q2_types.feature_data import (
    FeatureData, Sequence, Taxonomy, ProteinSequence
)
from q2_types.feature_table import FeatureTable, Frequency, PresenceAbsence
from q2_types.per_sample_sequences import (
    SequencesWithQuality, PairedEndSequencesWithQuality, MAGs, Contigs
)
from q2_types.sample_data import SampleData
from q2_types.feature_map import FeatureMap, MAGtoContigs
from qiime2.core.type import (
    Bool, Range, Int, Str, Float, List, Choices, Collection
)
from qiime2.core.type import (Properties, TypeMap)
from qiime2.plugin import (Plugin, Citations)
import q2_moshpit._examples as ex
import q2_moshpit
from q2_types.feature_data_mag import NOG, MAG
from q2_types.genome_data import (
    BLAST6, GenomeData, Loci, Genes, Proteins
)
from q2_types.kaiju import KaijuDB
from q2_types.kraken2 import (
    Kraken2Reports, Kraken2Outputs, Kraken2DB, Kraken2DBReport
)
from q2_types.kraken2._type import BrackenDB
from q2_types.per_sample_sequences._type import AlignmentMap
from q2_types.reference_db import (
    ReferenceDB, Diamond, Eggnog, NCBITaxonomy, EggnogProteinSequences
)

citations = Citations.load('citations.bib', package='q2_moshpit')

kraken2_params = {
    'threads': Int % Range(1, None),
    'confidence': Float % Range(0, 1, inclusive_end=True),
    'minimum_base_quality': Int % Range(0, None),
    'memory_mapping': Bool,
    'minimum_hit_groups': Int % Range(1, None),
    'quick': Bool,
    'report_minimizer_data': Bool
}
kraken2_param_descriptions = {
    'threads': 'Number of threads.',
    'confidence': 'Confidence score threshold.',
    'minimum_base_quality': 'Minimum base quality used in classification.'
                            ' Only applies when reads are used as input.',
    'memory_mapping': 'Avoids loading the database into RAM.',
    'minimum_hit_groups': 'Minimum number of hit groups (overlapping '
                          'k-mers sharing the same minimizer).',
    'quick': 'Quick operation (use first hit or hits).',
    'report_minimizer_data': 'Include number of read-minimizers per-taxon and'
                             ' unique read-minimizers per-taxon in the repot.'
}

partition_params = {"num_partitions": Int % Range(1, None)}
partition_param_descriptions = {
        "num_partitions": "The number of partitions to split the contigs"
        " into. Defaults to partitioning into individual"
        " samples."
}

plugin = Plugin(
    name='moshpit',
    version=q2_moshpit.__version__,
    website="https://github.com/bokulich-lab/q2-moshpit",
    package='q2_moshpit',
    description=(
        'MOdular SHotgun metagenome Pipelines with Integrated '
        'provenance Tracking: QIIME 2 plugin gor metagenome analysis with'
        'tools for genome binning and functional annotation.'),
    short_description='QIIME 2 plugin for metagenome analysis.',
)

importlib.import_module('q2_moshpit.eggnog')
importlib.import_module('q2_moshpit.metabat2')

plugin.methods.register_function(
    function=q2_moshpit.metabat2.bin_contigs_metabat,
    inputs={
        'contigs': SampleData[Contigs],
        'alignment_maps': SampleData[AlignmentMap]
    },
    parameters={
        'min_contig': Int % Range(1500, None),
        'max_p': Int % Range(1, 100),
        'min_s': Int % Range(1, 100),
        'max_edges': Int % Range(1, None),
        'p_tnf': Int % Range(0, 100),
        'no_add': Bool,
        'min_cv': Int % Range(1, None),
        'min_cv_sum': Int % Range(1, None),
        'min_cls_size': Int % Range(1, None),
        'num_threads': Int % Range(0, None),
        'seed': Int % Range(0, None),
        'debug': Bool,
        'verbose': Bool
    },
    outputs=[
        ('mags', SampleData[MAGs]),
        ('contig_map', FeatureMap[MAGtoContigs]),
        ('unbinned_contigs', SampleData[Contigs % Properties('unbinned')])
    ],
    input_descriptions={
        'contigs': 'Contigs to be binned.',
        'alignment_maps': 'Reads-to-contig alignment maps.'
    },
    parameter_descriptions={
        'min_contig': 'Minimum size of a contig for binning.',
        'max_p': 'Percentage of "good" contigs considered for binning '
                 'decided by connection among contigs. The greater, the '
                 'more sensitive.',
        'min_s': 'Minimum score of a edge for binning. The greater, the '
                 'more specific.',
        'max_edges': 'Maximum number of edges per node. The greater, the '
                     'more sensitive.',
        'p_tnf': 'TNF probability cutoff for building TNF graph. Use it to '
                 'skip the preparation step. (0: auto)',
        'no_add': 'Turning off additional binning for lost or small contigs.',
        'min_cv': 'Minimum mean coverage of a contig in each library '
                  'for binning.',
        'min_cv_sum': 'Minimum total effective mean coverage of a contig '
                      '(sum of depth over minCV) for binning.',
        'min_cls_size': 'Minimum size of a bin as the output.',
        'num_threads': 'Number of threads to use (0: use all cores).',
        'seed': 'For exact reproducibility. (0: use random seed)',
        'debug': 'Debug output.',
        'verbose': 'Verbose output.'
    },
    output_descriptions={
        'mags': 'The resulting MAGs.',
        'contig_map': 'Mapping of MAG identifiers to the contig identifiers '
                      'contained in each MAG.',
        'unbinned_contigs': 'Contigs that were not binned into any MAG.'
    },
    name='Bin contigs into MAGs using MetaBAT 2.',
    description='This method uses MetaBAT 2 to bin provided contigs '
                'into MAGs.',
    citations=[citations["kang2019"], citations["heng2009samtools"],
               citations["scikit_bio_release"]]
)

T_kraken_in, T_kraken_out_rep, T_kraken_out_hits = TypeMap({
    SampleData[SequencesWithQuality |
               PairedEndSequencesWithQuality]: (
        SampleData[Kraken2Reports % Properties('reads')],
        SampleData[Kraken2Outputs % Properties('reads')]
    ),
    SampleData[Contigs]: (
        SampleData[Kraken2Reports % Properties('contigs')],
        SampleData[Kraken2Outputs % Properties('contigs')]
    ),
    FeatureData[MAG]: (
        FeatureData[Kraken2Reports % Properties('mags')],
        FeatureData[Kraken2Outputs % Properties('mags')]
    ),
})

plugin.pipelines.register_function(
    function=q2_moshpit.kraken2.classification.classify_kraken2,
    inputs={
        "seqs": T_kraken_in,
        "kraken2_db": Kraken2DB,
    },
    parameters={**kraken2_params, **partition_params},
    outputs=[
        ('reports', T_kraken_out_rep),
        ('hits', T_kraken_out_hits),
    ],
    input_descriptions={
        "seqs": "Sequences to be classified. Single-/paired-end reads,"
                "contigs, or assembled MAGs can be provided.",
        "kraken2_db": "Kraken 2 database.",
    },
    parameter_descriptions={
        **kraken2_param_descriptions,
        **partition_param_descriptions
    },
    output_descriptions={
        'reports': 'Reports produced by Kraken2.',
        'hits': 'Output files produced by Kraken2.',
    },
    name='Perform taxonomic classification of reads or MAGs using Kraken 2.',
    description='Use Kraken 2 to classify provided DNA sequence reads, '
                'contigs, or MAGs into taxonomic groups.',
    citations=[citations["wood2019"]]
)

plugin.methods.register_function(
    function=q2_moshpit.kraken2.classification._classify_kraken2,
    inputs={
        "seqs": T_kraken_in,
        "kraken2_db": Kraken2DB,
    },
    parameters=kraken2_params,
    outputs=[
        ('reports', T_kraken_out_rep),
        ('hits', T_kraken_out_hits),
    ],
    input_descriptions={
        "seqs": "Sequences to be classified. Single-/paired-end reads,"
                "contigs, or assembled MAGs can be provided.",
        "kraken2_db": "Kraken 2 database.",
    },
    parameter_descriptions=kraken2_param_descriptions,
    output_descriptions={
        'reports': 'Reports produced by Kraken2.',
        'hits': 'Output files produced by Kraken2.',
    },
    name='Perform taxonomic classification of reads or MAGs using Kraken 2.',
    description='Use Kraken 2 to classify provided DNA sequence reads, '
                'contigs, or MAGs into taxonomic groups.',
    citations=[citations["wood2019"]]
)

T_kraken_collate_reports_in, T_kraken_collate_reports_out = TypeMap({
    SampleData[Kraken2Reports % Properties('reads', 'contigs')]: (
        SampleData[Kraken2Reports % Properties('reads', 'contigs')],
    ),
    SampleData[Kraken2Reports % Properties('reads')]: (
        SampleData[Kraken2Reports % Properties('reads')],
    ),
    SampleData[Kraken2Reports % Properties('contigs')]: (
        SampleData[Kraken2Reports % Properties('contigs')],
    )
})

plugin.methods.register_function(
    function=q2_moshpit.kraken2.helpers.collate_kraken2_reports,
    inputs={"kraken2_reports": List[T_kraken_collate_reports_in]},
    parameters={},
    outputs={"collated_kraken2_reports": T_kraken_collate_reports_out},
    name="Collate kraken2 reports",
    description="Collates kraken2 reports"
)

T_kraken_collate_outputs_in, T_kraken_collate_outputs_out = TypeMap({
    SampleData[Kraken2Outputs % Properties('reads', 'contigs')]: (
        SampleData[Kraken2Outputs % Properties('reads', 'contigs')],
    ),
    SampleData[Kraken2Outputs % Properties('reads')]: (
        SampleData[Kraken2Outputs % Properties('reads')],
    ),
    SampleData[Kraken2Outputs % Properties('contigs')]: (
        SampleData[Kraken2Outputs % Properties('contigs')],
    )
})

plugin.methods.register_function(
    function=q2_moshpit.kraken2.helpers.collate_kraken2_outputs,
    inputs={"kraken2_outputs": List[T_kraken_collate_outputs_in]},
    parameters={},
    outputs={"collated_kraken2_outputs": T_kraken_collate_outputs_out},
    name="Collate kraken2 outputs.",
    description="Collates kraken2 outputs."
)

plugin.methods.register_function(
    function=q2_moshpit.kraken2.bracken.estimate_bracken,
    inputs={
        "kraken_reports": SampleData[Kraken2Reports % Properties('reads')],
        "bracken_db": BrackenDB
    },
    parameters={
        'threshold': Int % Range(0, None),
        'read_len': Int % Range(0, None),
        'level': Str % Choices(['D', 'P', 'C', 'O', 'F', 'G', 'S'])
    },
    outputs=[
        ('reports', SampleData[Kraken2Reports % Properties('bracken')]),
        ('taxonomy', FeatureData[Taxonomy]),
        ('table', FeatureTable[Frequency])
    ],
    input_descriptions={
        "kraken_reports": "Reports produced by Kraken2.",
        "bracken_db": "Bracken database."
    },
    parameter_descriptions={
        'threshold': 'Bracken: number of reads required PRIOR to abundance '
                     'estimation to perform re-estimation.',
        'read_len': 'Bracken: the ideal length of reads in your sample. '
                    'For paired end data (e.g., 2x150) this should be set '
                    'to the length of the single-end reads (e.g., 150).',
        'level': 'Bracken: specifies the taxonomic rank to analyze. Each '
                 'classification at this specified rank will receive an '
                 'estimated number of reads belonging to that rank after '
                 'abundance estimation.'
    },
    output_descriptions={
        'reports': 'Reports modified by Bracken.',
    },
    name='Perform read abundance re-estimation using Bracken.',
    description='This method uses Bracken to re-estimate read abundances.',
    citations=[citations["wood2019"]]
)

plugin.methods.register_function(
    function=q2_moshpit.kraken2.build_kraken_db,
    inputs={
        "seqs": List[FeatureData[Sequence]]
    },
    parameters={
        'collection': Str % Choices(
            ['viral', 'minusb', 'standard', 'standard8',
             'standard16', 'pluspf', 'pluspf8', 'pluspf16',
             'pluspfp', 'pluspfp8', 'pluspfp16', 'eupathdb',
             'nt', 'greengenes', 'rdp', 'silva132', 'silva138'],
        ),
        'threads': Int % Range(1, None),
        'kmer_len': Int % Range(1, None),
        'minimizer_len': Int % Range(1, None),
        'minimizer_spaces': Int % Range(1, None),
        'no_masking': Bool,
        'max_db_size': Int % Range(0, None),
        'use_ftp': Bool,
        'load_factor': Float % Range(0, 1),
        'fast_build': Bool,
        'read_len': List[Int % Range(1, None)],
    },
    outputs=[
        ('kraken2_database', Kraken2DB),
        ('bracken_database', BrackenDB),
    ],
    input_descriptions={
        "seqs": "Sequences to be added to the Kraken 2 database."
    },
    parameter_descriptions={
        'collection': 'Name of the database collection to be fetched. '
                      'Please check https://benlangmead.github.io/aws-'
                      'indexes/k2 for the description of the available '
                      'options.',
        'threads': 'Number of threads. Only applicable when building a '
                   'custom database.',
        'kmer_len': 'K-mer length in bp/aa.',
        'minimizer_len': 'Minimizer length in bp/aa.',
        'minimizer_spaces': 'Number of characters in minimizer that are '
                            'ignored in comparisons.',
        'no_masking': 'Avoid masking low-complexity sequences prior to '
                      'building; masking requires dustmasker or segmasker '
                      'to be installed in PATH.',
        'max_db_size': 'Maximum number of bytes for Kraken 2 hash table; '
                       'if the estimator determines more would normally be '
                       'needed, the reference library will be downsampled '
                       'to fit.',
        'use_ftp': 'Use FTP for downloading instead of RSYNC.',
        'load_factor': 'Proportion of the hash table to be populated.',
        'fast_build': 'Do not require database to be deterministically '
                      'built when using multiple threads. This is faster, '
                      'but does introduce variability in minimizer/LCA pairs.',
        'read_len': 'Ideal read lengths to be used while building the Bracken '
                    'database.'
    },
    output_descriptions={
        'kraken2_database': 'Kraken2 database.',
        'bracken_database': 'Bracken database.'
    },
    name='Build Kraken 2 database.',
    description='This method builds Kraken 2 and Bracken databases either (1) '
                'from provided DNA sequences to build a custom database, or '
                '(2) simply fetches pre-built versions from an online '
                'resource.',
    citations=[citations["wood2019"], citations["lu2017"]]
)

plugin.methods.register_function(
    function=q2_moshpit.kraken2.database.inspect_kraken2_db,
    inputs={"db": Kraken2DB},
    parameters={"threads": Int % Range(1, None)},
    outputs=[("report", Kraken2DBReport)],
    input_descriptions={
        "db": "The Kraken 2 database for which to generate the report."
    },
    parameter_descriptions={
        "threads": "The number of threads to use."
    },
    output_descriptions={
        "report": "The report of the supplied database."
    },
    name="Inspect a Kraken 2 database.",
    description="This method generates a report of identical format to those "
                "generated by classify_kraken2, with a slightly different "
                "interpretation. Instead of reporting the number of inputs "
                "classified to a taxon/clade, the report displays the number "
                "of minimizers mapped to each taxon/clade.",
    citations=[citations["wood2019"]],
)

plugin.methods.register_function(
    function=q2_moshpit.dereplication.dereplicate_mags,
    inputs={
        "mags": SampleData[MAGs],
        "distance_matrix": DistanceMatrix
    },
    parameters={
        "threshold": Float % Range(0, 1, inclusive_end=True)
    },
    outputs=[
        ('dereplicated_mags', FeatureData[MAG]),
        ('feature_table', FeatureTable[PresenceAbsence])
    ],
    input_descriptions={
        "mags": "MAGs to be dereplicated.",
        "distance_matrix": "Matrix of distances between MAGs."
    },
    parameter_descriptions={
        "threshold": "Similarity threshold required to consider "
                     "two bins identical."
    },
    output_descriptions={
        "dereplicated_mags": "Dereplicated MAGs.",
        "feature_table": "Mapping between MAGs and samples."
    },
    name='Dereplicate MAGs from multiple samples.',
    description='This method dereplicates MAGs from multiple samples '
                'using distances between them found in the provided '
                'distance matrix. For each cluster of similar MAGs, '
                'the longest one will be selected as the representative.',
    citations=[]
)

select_features_taxonomy_description = (
    'Output taxonomy. Infra-clade ranks are ignored unless if they are '
    'strain-level. Missing internal ranks are annotated by their next '
    'most specific rank, with the exception of k__Bacteria and k__Archaea, '
    'which match their domain name.')

select_features_description = (
    'Convert a Kraken 2 report, which is an annotated NCBI taxonomy tree, '
    'into generic artifacts for downstream analyses.')

plugin.methods.register_function(
    function=q2_moshpit.kraken2.kraken2_to_features,
    inputs={
        'reports': SampleData[Kraken2Reports]
    },
    parameters={
        'coverage_threshold': Float % Range(0, 100, inclusive_end=True)
    },
    outputs=[
        ('table', FeatureTable[PresenceAbsence]),
        ('taxonomy', FeatureData[Taxonomy])
    ],
    input_descriptions={
        'reports': 'Per-sample Kraken 2 reports.'
    },
    parameter_descriptions={
        'coverage_threshold': 'The minimum percent coverage required to '
                              'produce a feature.'
    },
    output_descriptions={
        'table': 'A presence/absence table of selected features. The features '
                 'are not of even ranks, but will be the most specific rank '
                 'available.',
        'taxonomy': select_features_taxonomy_description,
    },
    name='Select features from a Kraken 2 report.',
    description=select_features_description
)

plugin.methods.register_function(
    function=q2_moshpit.kraken2.kraken2_to_mag_features,
    inputs={
        'reports': FeatureData[Kraken2Reports % Properties('mags')],
        'hits': FeatureData[Kraken2Outputs % Properties('mags')],
    },
    parameters={
        'coverage_threshold': Float % Range(0, 100, inclusive_end=True),
        # 'lca_mode': Str % Choices(['lca', 'majority'])
    },
    outputs=[('taxonomy', FeatureData[Taxonomy])],
    input_descriptions={
        'reports': 'Per-sample Kraken 2 reports.',
        'hits': 'Per-sample Kraken 2 output files.',
    },
    parameter_descriptions={
        'coverage_threshold': 'The minimum percent coverage required to '
                              'produce a feature.',
        # 'lca_mode': 'The method used to determine the LCA of a MAG using '
        #             'taxonomic assignments of its contigs. '
    },
    output_descriptions={
        'taxonomy': select_features_taxonomy_description,
    },
    name='Select MAG features from a Kraken 2 report.',
    description=select_features_description
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog.build_custom_diamond_db,
    inputs={
        'seqs': FeatureData[ProteinSequence],
        'taxonomy': ReferenceDB[NCBITaxonomy],
    },
    input_descriptions={
        'seqs': "Protein reference database.",
        'taxonomy': "Reference taxonomy, "
                    "needed to provide taxonomy features."
    },
    outputs=[('diamond_db', ReferenceDB[Diamond])],
    output_descriptions={
        'diamond_db': "DIAMOND database."
    },
    parameters={
        "threads": Int % Range(1, None),
        "file_buffer_size": Int % Range(1, None),
        "ignore_warnings": Bool,
        "no_parse_seqids": Bool
    },
    parameter_descriptions={
        "threads": "Number of CPU threads.",
        "file_buffer_size": "File buffer size in bytes.",
        "ignore_warnings": "Ignore warnings.",
        "no_parse_seqids": "Print raw seqids without parsing."
    },
    name="Create a DIAMOND formatted reference database from a FASTA input "
         "file.",
    description="Creates an artifact containing a binary DIAMOND database "
                "file (ref_db.dmnd) from a protein reference database "
                "file in FASTA format.",
    citations=[citations["buchfink_sensitive_2021"]],
    examples={
        "Minimum working example": ex.diamond_makedb
    }
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog.fetch_eggnog_db,
    inputs={},
    parameters={},
    outputs=[("eggnog_db", ReferenceDB[Eggnog])],
    output_descriptions={
        "eggnog_db": "eggNOG annotation database."
    },
    name="Fetch the databases necessary to run the "
         "eggnog-annotate action.",
    description="Downloads eggnog reference database  "
                "using the `download_eggnog_data.py` script from eggNOG. "
                "Here, this script downloads 3 files "
                "and creates an artifact with them. At least 80 Gb of "
                "storage space is required to run this action. "
                "Links to files: "
                "eggnog.db: "
                "http://eggnogdb.embl.de/download/emapperdb-5.0.2/eggnog.db.gz"
                "eggnog.taxa.db: "
                "http://eggnogdb.embl.de/download/emapperdb-5.0.2/"
                "eggnog.taxa.tar.gz"
                "eggnog.taxa.db.traverse.pkl: "
                "http://eggnogdb.embl.de/download/emapperdb-5.0.2/"
                "eggnog_proteins.dmnd.gz",
    citations=[citations["huerta_cepas_eggnog_2019"]]
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog.fetch_diamond_db,
    inputs={},
    parameters={},
    outputs=[("diamond_db", ReferenceDB[Diamond])],
    output_descriptions={
        "diamond_db": "Complete Diamond reference database."
    },
    name="Fetch the complete Diamond database necessary to run the "
         "eggnog-diamond-search action.",
    description="Downloads Diamond reference database.  "
                "This action downloads 1 file (ref_db.dmnd). "
                "At least 18 Gb of storage space is "
                "required to run this action. "
                "Link to database: "
                "http://eggnogdb.embl.de/download/emapperdb-5.0.2/"
                "eggnog_proteins.dmnd.gz",
    citations=[
        citations["buchfink_sensitive_2021"],
        citations["huerta_cepas_eggnog_2019"]
    ]
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog.fetch_eggnog_proteins,
    inputs={},
    parameters={},
    outputs=[("eggnog_proteins", ReferenceDB[EggnogProteinSequences])],
    output_descriptions={
        "eggnog_proteins": "eggNOG database of protein sequences and "
                           "their corresponding taxonomy information."
    },
    name="Fetch the databases necessary to run the "
         "build-eggnog-diamond-db action.",
    description="Downloads eggnog proteome database.  "
                "This script downloads 2 files "
                "(e5.proteomes.faa and e5.taxid_info.tsv) "
                "and creates and artifact with them. At least 18 GB of "
                "storage space is required to run this action. ",
    citations=[citations["huerta_cepas_eggnog_2019"]]
)


plugin.methods.register_function(
    function=q2_moshpit.eggnog.fetch_ncbi_taxonomy,
    inputs={},
    parameters={},
    outputs=[("taxonomy", ReferenceDB[NCBITaxonomy])],
    output_descriptions={
        "taxonomy": "NCBI reference taxonomy."
    },
    name="Fetch NCBI reference taxonomy.",
    description="Downloads NCBI reference taxonomy from the NCBI FTP server. "
                "The resulting artifact is required by the "
                "build-custom-diamond-db action if one wishes to "
                "create a Diamond data base with taxonomy features. "
                "At least 30 GB of "
                "storage space is required to run this action.",
    citations=[citations["NCBI"]]
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog.build_eggnog_diamond_db,
    inputs={
        'eggnog_proteins': ReferenceDB[EggnogProteinSequences],
    },
    input_descriptions={
        'eggnog_proteins': "eggNOG database of protein sequences and "
                           "their corresponding taxonomy information "
                           "(generated through the `fetch-eggnog-proteins` "
                           "action)."
    },
    parameters={
        'taxon': Int % Range(2, 1579337)
    },
    parameter_descriptions={
        'taxon': "NCBI Taxon ID number."
    },
    outputs=[("diamond_db", ReferenceDB[Diamond])],
    output_descriptions={
        "diamond_db": "Complete Diamond reference database for the "
                      "specified taxon."
    },
    name="Create a DIAMOND formatted reference database for the "
         "specified taxon.",
    description="Creates a DIAMOND database which contains the protein "
                "sequences that belong to the specified taxon.",
    citations=[
        citations["buchfink_sensitive_2021"],
        citations["huerta_cepas_eggnog_2019"]
    ]
)

plugin.pipelines.register_function(
    function=q2_moshpit.eggnog.eggnog_diamond_search,
    inputs={
        'sequences': SampleData[Contigs] | FeatureData[MAG],
        'diamond_db': ReferenceDB[Diamond],
    },
    parameters={
        'num_cpus': Int,
        'db_in_memory': Bool,
        **partition_params
    },
    input_descriptions={
        'sequences': 'Sequence data of the contigs we want to '
                     'search for hits using the Diamond Database',
        'diamond_db': 'The filepath to an artifact containing the '
                      'Diamond database',
    },
    parameter_descriptions={
        'num_cpus': 'Number of CPUs to utilize. \'0\' will '
                    'use all available.',
        'db_in_memory': 'Read database into memory. The '
                        'database can be very large, so this '
                        'option should only be used on clusters or other '
                        'machines with enough memory.',
        **partition_param_descriptions
    },
    outputs=[
        ('eggnog_hits', SampleData[BLAST6]),
        ('table', FeatureTable[Frequency])
    ],
    name='Run eggNOG search using diamond aligner',
    description="This method performs the steps by which we find our "
                "possible target sequences to annotate using the diamond "
                "search functionality from the eggnog `emapper.py` script",
    citations=[
        citations["buchfink_sensitive_2021"],
        citations["huerta_cepas_eggnog_2019"]
    ]
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog._eggnog_diamond_search,
    inputs={
        'sequences': SampleData[Contigs] | FeatureData[MAG],
        'diamond_db': ReferenceDB[Diamond],
    },
    parameters={
        'num_cpus': Int,
        'db_in_memory': Bool,
    },
    input_descriptions={
        'sequences': 'Contigs or MAGs to '
                     'search against the Diamond Database.',
        'diamond_db': 'The Diamond database.',
    },
    parameter_descriptions={
        'num_cpus': 'Number of CPUs to utilize. \'0\' will '
                    'use all available CPUs.',
        'db_in_memory': 'Read database into memory. The '
                        'database can be very large, so this '
                        'option should only be used on clusters or other '
                        'machines with enough memory.',
    },
    outputs=[
        ('eggnog_hits', SampleData[BLAST6]),
        ('table', FeatureTable[Frequency])
    ],
    name='Run eggNOG search using diamond aligner.',
    description="Use Diamond and eggNOG to align contig or MAG sequences "
                "against the Diamond database.",
    citations=[
        citations["buchfink_sensitive_2021"],
        citations["huerta_cepas_eggnog_2019"]
    ]
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog._eggnog_feature_table,
    inputs={
        'seed_orthologs': SampleData[BLAST6]
    },
    parameters={},
    input_descriptions={
        'seed_orthologs': 'Sequence data to be turned into an eggnog feature '
                          'table.'
    },
    parameter_descriptions={},
    outputs=[
        ('table', FeatureTable[Frequency])
    ],
    name='Create an eggnog table',
    description='Create an eggnog table'
)

plugin.pipelines.register_function(
    function=q2_moshpit.eggnog.eggnog_annotate,
    inputs={
        'eggnog_hits': SampleData[BLAST6],
        'eggnog_db': ReferenceDB[Eggnog],
    },
    parameters={
        'db_in_memory': Bool,
        'num_cpus': Int % Range(0, None),
        **partition_params
    },
    parameter_descriptions={
        'db_in_memory': 'Read eggnog database into memory. The '
                        'eggnog database is very large(>44GB), so this '
                        'option should only be used on clusters or other '
                        'machines with enough memory.',
        'num_cpus': 'Number of CPUs to utilize. \'0\' will '
                    'use all available.',
        **partition_param_descriptions
    },
    outputs=[('ortholog_annotations', FeatureData[NOG])],
    name='Annotate orthologs against eggNOG database.',
    description="Apply eggnog mapper to annotate seed orthologs.",
    citations=[citations["huerta_cepas_eggnog_2019"]]
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog._eggnog_annotate,
    inputs={
        'eggnog_hits': SampleData[BLAST6],
        'eggnog_db': ReferenceDB[Eggnog],
    },
    parameters={
        'db_in_memory': Bool,
        'num_cpus': Int % Range(0, None)
    },
    parameter_descriptions={
        'db_in_memory': 'Read eggnog database into memory. The '
                        'eggnog database is very large(>44GB), so this '
                        'option should only be used on clusters or other '
                        'machines with enough memory.',
        'num_cpus': 'Number of CPUs to utilize. \'0\' will '
                    'use all available.',
    },
    outputs=[('ortholog_annotations', FeatureData[NOG])],
    name='Annotate orthologs against eggNOG database',
    description="Apply eggnog mapper to annotate seed orthologs.",
    citations=[citations["huerta_cepas_eggnog_2019"]]
)

plugin.methods.register_function(
    function=q2_moshpit.partition.partition_sample_data_mags,
    inputs={"mags": SampleData[MAGs]},
    parameters={"num_partitions": Int % Range(1, None)},
    outputs={"partitioned_mags": Collection[SampleData[MAGs]]},
    input_descriptions={"mags": "The MAGs to partition."},
    parameter_descriptions={
        "num_partitions": "The number of partitions to split the MAGs"
        " into. Defaults to partitioning into individual"
        " MAGs."
    },
    name="Partition MAGs",
    description="Partition a SampleData[MAGs] artifact into smaller "
                "artifacts containing subsets of the MAGs",
)

plugin.methods.register_function(
    function=q2_moshpit.partition.partition_orthologs,
    inputs={"orthologs": SampleData[BLAST6]},
    parameters={"num_partitions": Int % Range(1, None)},
    outputs={"partitioned_orthologs": Collection[SampleData[BLAST6]]},
    input_descriptions={"orthologs": "The orthologs to partition."},
    parameter_descriptions={
        "num_partitions": "The number of partitions to split the MAGs"
        " into. Defaults to partitioning into individual"
        " MAGs."
    },
    name="Partition orthologs",
    description="Partition a SampleData[MAGs] artifact into smaller "
                "artifacts containing subsets of the MAGs",
)

plugin.methods.register_function(
    function=q2_moshpit.partition.collate_sample_data_mags,
    inputs={"mags": List[SampleData[MAGs]]},
    parameters={},
    outputs={"collated_mags": SampleData[MAGs]},
    input_descriptions={"mags": "A collection of MAGs to be collated."},
    name="Collate mags",
    description="Takes a collection of SampleData[MAGs]'s "
                "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=q2_moshpit.partition.partition_feature_data_mags,
    inputs={"mags": FeatureData[MAG]},
    parameters={"num_partitions": Int % Range(1, None)},
    outputs={"partitioned_mags": Collection[FeatureData[MAG]]},
    input_descriptions={"mags": "MAGs to partition."},
    parameter_descriptions={
        "num_partitions": "The number of partitions to split the MAGs"
        " into. Defaults to partitioning into individual"
        " MAGs."
    },
    name="Partition MAGs",
    description="Partition a FeatureData[MAG] artifact into smaller "
                "artifacts containing subsets of the MAGs",
)

plugin.methods.register_function(
    function=q2_moshpit.partition.collate_feature_data_mags,
    inputs={"mags": List[FeatureData[MAG]]},
    parameters={},
    outputs={"collated_mags": FeatureData[MAG]},
    input_descriptions={"mags": "A collection of MAGs to be collated."},
    name="Collate mags",
    description="Takes a collection of FeatureData[MAG]'s "
                "and collates them into a single artifact.",
)

plugin.methods.register_function(
    function=q2_moshpit.partition.collate_orthologs,
    inputs={"orthologs": List[SampleData[BLAST6]]},
    parameters={},
    outputs={"collated_orthologs": SampleData[BLAST6]},
    input_descriptions={"orthologs": "Orthologs to collate"},
    parameter_descriptions={},
    name="Collate Orthologs",
    description="Collate a List of SampleData[BLAST6] into one"
)

plugin.methods.register_function(
    function=q2_moshpit.partition.collate_annotations,
    inputs={'ortholog_annotations': List[FeatureData[NOG]]},
    parameters={},
    outputs=[('collated_annotations', FeatureData[NOG])],
    input_descriptions={
        'ortholog_annotations': "Collection of ortholog annotations."
    },
    output_descriptions={
        'collated_annotations': "Collated ortholog annotations."
    },
    name='Collate ortholog annotations.',
    description="Takes a collection of FeatureData[NOG]'s "
                "and collates them into a single artifact.",
)

busco_params = {
    "mode": Str % Choices(["genome"]),
    "lineage_dataset": Str,
    "augustus": Bool,
    "augustus_parameters": Str,
    "augustus_species": Str,
    "auto_lineage": Bool,
    "auto_lineage_euk": Bool,
    "auto_lineage_prok": Bool,
    "cpu": Int % Range(1, None),
    "config": Str,
    "contig_break": Int % Range(0, None),
    "evalue": Float % Range(0, None, inclusive_start=False),
    "force": Bool,
    "limit": Int % Range(1, 20),
    "long": Bool,
    "metaeuk_parameters": Str,
    "metaeuk_rerun_parameters": Str,
    "miniprot": Bool,
    "scaffold_composition": Bool,
}
busco_param_descriptions = {
    "mode": "Specify which BUSCO analysis mode to run."
            "Currently only the 'genome' option is supported, "
            "for genome assemblies. In the future modes for transcriptome "
            "assemblies and for annotated gene sets (proteins) will be made "
            "available.",
    "lineage_dataset": "Specify the name of the BUSCO lineage to be used. "
                       "To see all possible options run `busco "
                       "--list-datasets`.",
    "augustus": "Use augustus gene predictor for eukaryote runs.",
    "augustus_parameters": "Pass additional arguments to Augustus. "
                           "All arguments should be contained within a single "
                           "string with no white space, with each argument "
                           "separated by a comma. "
                           "Example: '--PARAM1=VALUE1,--PARAM2=VALUE2'.",
    "augustus_species": "Specify a species for Augustus training.",
    "auto_lineage": "Run auto-lineage to find optimum lineage path.",
    "auto_lineage_euk": "Run auto-placement just on eukaryote tree to find "
                        "optimum lineage path.",
    "auto_lineage_prok": "Run auto-lineage just on non-eukaryote trees to "
                         "find optimum lineage path.",
    "cpu": "Specify the number (N=integer) of threads/cores to use.",
    "config": "Provide a config file.",
    "contig_break": "Number of contiguous Ns to signify a break between "
                    "contigs. Default is n=10. "
                    "See https://gitlab.com/ezlab/busco/-/issues/691 for a "
                    "more detailed explanation.",
    "evalue": "E-value cutoff for BLAST searches. "
              "Allowed formats, 0.001 or 1e-03, Default: 1e-03.",
    "force": "Force rewriting of existing files. Must be used when output "
             "files with the provided name already exist.",
    "limit": "How many candidate regions (contig or transcript) to consider "
             "per BUSCO. Default: 3.",
    "long": "Optimization Augustus self-training mode (Default: Off); "
            "adds considerably to the run time, "
            "but can improve results for some non-model organisms.",
    "metaeuk_parameters": "Pass additional arguments to Metaeuk for the first "
                          "run. All arguments should be contained within a "
                          "single string with no white space, with each "
                          "argument separated by a comma. "
                          "Example: `--PARAM1=VALUE1,--PARAM2=VALUE2`.",
    "metaeuk_rerun_parameters": "Pass additional arguments to Metaeuk for the "
                                "second run. All arguments should be "
                                "contained within a single string with no "
                                "white space, with each argument separated by "
                                "a comma. "
                                "Example: `--PARAM1=VALUE1,--PARAM2=VALUE2`.",
    "miniprot": "Use miniprot gene predictor for eukaryote runs.",
    "scaffold_composition": "Writes ACGTN content per scaffold to a file "
                            "`scaffold_composition.txt`.",
}


plugin.visualizers.register_function(
    function=q2_moshpit.busco.evaluate_busco,
    inputs={
        "bins": SampleData[MAGs],
    },
    parameters=busco_params,
    input_descriptions={
        "bins": "MAGs to be analyzed.",
    },
    parameter_descriptions=busco_param_descriptions,
    name="Evaluate quality of the generated MAGs using BUSCO.",
    description="This method uses BUSCO "
                "(Benchmarking Universal Single-Copy Ortholog assessment "
                "tool) to assess the quality of assembled MAGs and generates "
                "visualizations summarizing the results.",
    citations=[citations["manni_busco_2021"]],
)

plugin.methods.register_function(
    function=q2_moshpit.prodigal.predict_genes_prodigal,
    inputs={
        'mags': FeatureData[MAG]
    },
    input_descriptions={
        'mags': 'MAGs for which one wishes to predict genes.'
    },
    parameters={
        "translation_table_number": Str % Choices([
            '1', '2', '3', '4', '5', '6',
            '9', '10', '11', '12', '13', '14', '15', '16',
            '21', '22', '23', '24', '25'
        ])
    },
    parameter_descriptions={
        'translation_table_number':
            'Translation table to be used to '
            'translate genes into a sequence of amino '
            'acids. See '
            'https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi '
            'for reference.'
    },
    outputs=[
        ('loci', GenomeData[Loci]),
        ('genes', GenomeData[Genes]),
        ('proteins', GenomeData[Proteins])
    ],
    output_descriptions={
        'loci': "Gene coordinates files (one per MAG) listing the location of "
                "each predicted gene as well as some additional scoring "
                "information.",
        'genes': "Fasta files (one per MAG) with the nucleotide sequences of "
                 "the predicted genes.",
        'proteins': "Fasta files (one per MAG) with the protein translation "
                    "of the predicted genes."
    },
    name='Predict gene sequences from MAGs using Prodigal.',
    description="Prodigal (PROkaryotic DYnamic programming "
                "Gene-finding ALgorithm), a gene prediction algorithm "
                "designed for improved gene structure prediction, translation "
                "initiation site recognition, and reduced false positives in "
                "bacterial and archaeal genomes.",
    citations=[citations["hyatt_prodigal_2010"]]
)

plugin.methods.register_function(
    function=q2_moshpit.kaiju.fetch_kaiju_db,
    inputs={},
    parameters={
        "database_type": Str
        % Choices(
            [
                "nr",
                "nr_euk",
                "refseq",
                "fungi",
                "viruses",
                "plasmids",
                "progenomes",
                "rvdb",
            ]
        ),
    },
    outputs=[
        ("database", KaijuDB),
    ],
    input_descriptions={},
    parameter_descriptions={
        "database_type": "Type of database to be downloaded. For more "
        "information on available types please see the list on "
        "Kaiju's web server: https://kaiju.binf.ku.dk/server",
    },
    output_descriptions={"database": "Kaiju database."},
    name="Fetch Kaiju database.",
    description="This method fetches the latest Kaiju database from "
                "https://kaiju.binf.ku.dk/server.",
    citations=[citations["menzel2016"]],
)

plugin.methods.register_function(
    function=q2_moshpit.kaiju.classify_kaiju,
    inputs={
        "seqs": SampleData[
            SequencesWithQuality | PairedEndSequencesWithQuality
            ],
        "db": KaijuDB,
    },
    parameters={
        "z": Int % Range(1, None),
        "a": Str % Choices(["greedy", "mem"]),
        "e": Int % Range(1, None),
        "m": Int % Range(1, None),
        "s": Int % Range(1, None),
        "evalue": Float % Range(0, 1),
        "x": Bool,
        "r": Str % Choices(
            ["phylum", "class", "order", "family", "genus", "species"]
        ),
        "c": Float % Range(0, 100, inclusive_start=True),
        "exp": Bool,
        "u": Bool,
    },
    outputs=[
        ("abundances", FeatureTable[Frequency]),
        ("taxonomy", FeatureData[Taxonomy])
    ],
    input_descriptions={
        "seqs": "Sequences to be classified.",
        "db": "Kaiju database.",
    },
    parameter_descriptions={
        "z": "Number of threads.",
        "a": "Run mode.",
        "e": "Number of mismatches allowed in Greedy mode.",
        "m": "Minimum match length.",
        "s": "Minimum match score in Greedy mode.",
        "evalue": "Minimum E-value in Greedy mode.",
        "x": "Enable SEG low complexity filter.",
        "r": "Taxonomic rank.",
        "c": "Minimum required number or fraction of reads for "
             "the taxon  (except viruses) to be reported.",
        "exp": "Expand viruses, which are always shown as full "
               "taxon path and read counts are not summarized in "
               "higher taxonomic levels.",
        "u": "Do not count unclassified reads for the total reads "
             "when calculating percentages for classified reads."
    },
    output_descriptions={
        "abundances": "Read abundances.",
        "taxonomy": "Linked taxonomy."
    },
    name="Classify reads using Kaiju.",
    description="This method uses Kaiju to perform taxonomic "
                "classification of DNA sequence reads.",
    citations=[citations["menzel2016"]],
)
