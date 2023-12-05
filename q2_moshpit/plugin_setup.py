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
    SequencesWithQuality, PairedEndSequencesWithQuality
)
from q2_types.sample_data import SampleData
from q2_types.feature_map import FeatureMap, MAGtoContigs
from qiime2.core.type import Bool, Range, Int, Str, Float, List, Choices
from qiime2.core.type import (Properties, TypeMap)
from qiime2.plugin import (Plugin, Citations)
import q2_moshpit._examples as ex
import q2_moshpit
from q2_types_genomics.feature_data import NOG, MAG
from q2_types_genomics.genome_data import (
    BLAST6, GenomeData, Loci, Genes, Proteins
)
from q2_types_genomics.kaiju import KaijuDB
from q2_types_genomics.kraken2 import (
    Kraken2Reports, Kraken2Outputs, Kraken2DB, Kraken2DBReport
)
from q2_types_genomics.kraken2._type import BrackenDB
from q2_types_genomics.per_sample_data import MAGs, Contigs
from q2_types_genomics.per_sample_data._type import AlignmentMap
from q2_types_genomics.reference_db import (
    ReferenceDB, Diamond, Eggnog, TaxonomyNCBI
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
        'contigs': 'Placeholder.',
        'alignment_maps': 'Placeholder.'
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
    citations=[citations["kang2019"]]
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

plugin.methods.register_function(
    function=q2_moshpit.kraken2.classification.classify_kraken2,
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
        "seqs": "The sequences to be classified. Single-end or paired-end "
                "reads, contigs, or MAGs can be provided.",
        "kraken2_db": "Kraken 2 database.",
    },
    parameter_descriptions=kraken2_param_descriptions,
    output_descriptions={
        'reports': 'Reports produced by Kraken2.',
        'hits': 'Output files produced by Kraken2.',
    },
    name='Perform taxonomic classification of reads or MAGs using Kraken 2.',
    description='This method uses Kraken 2 to classify provided NGS reads '
                'or MAGs into taxonomic groups.',
    citations=[citations["wood2019"]]
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
        'read_len': 'Bracken: read length to get all classifications for.',
        'level': 'Bracken: taxonomic level to estimate abundance at.'
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
             'pluspfp', 'pluspfp8', 'pluspfp16', 'eupathdb'],
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
                      'to be installed in PATH',
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
    description='This method builds a Kraken 2/Bracken databases from '
                'provided DNA sequences or simply fetches pre-built '
                'versions from an online resource.',
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
                "generated by classify_kraken2, with the interpretation being "
                "instead of reporting the number of inputs classified to a "
                "taxon/clade, the number of minimizers mapped to a "
                "taxon/clade are reported.",
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
        'coverage_threshold': 'The minimum percent coverage required to'
                              ' produce a feature.'
    },
    output_descriptions={
        'table': 'A presence/absence table of selected features. The features'
                 ' are not of even ranks, but will be the most specific rank'
                 ' available.',
        'taxonomy': 'Infra-clade ranks are ignored '
                    'unless they are strain-level. Missing internal ranks '
                    'are annotated by their next most specific rank, '
                    'with the exception of k__Bacteria and k__Archaea which '
                    'match their domain\'s name.',
    },
    name='Select downstream features from Kraken 2',
    description='Convert a Kraken 2 report, which is an annotated NCBI '
                'taxonomy tree into generic artifacts for downstream '
                'analyses.'
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
        'taxonomy': 'Infra-clade ranks are ignored '
                    'unless they are strain-level. Missing internal ranks '
                    'are annotated by their next most specific rank, '
                    'with the exception of k__Bacteria and k__Archaea which '
                    'match their domain\'s name.',
    },
    name='Select downstream MAG features from Kraken 2',
    description='Convert a Kraken 2 report, which is an annotated NCBI '
                'taxonomy tree into generic artifacts for downstream '
                'analyses.'
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog.build_custom_diamond_db,
    inputs={
        'sequences': FeatureData[ProteinSequence],
        'taxonomy_data': ReferenceDB[TaxonomyNCBI],
    },
    input_descriptions={
        'sequences': "Artifact containing protein reference database file "
                     "in FASTA format.",
        'taxonomy_data': "Artifact containing taxonomy data. "
                         "Needed in order to provide taxonomy features. "
                         "Can be generated through name_of_action."
                         # TODO: update action name here
    },
    outputs=[('diamond_db', ReferenceDB[Diamond])],
    output_descriptions={
        'diamond_db': "Artifact containing a binary DIAMOND database file."
    },
    parameters={
        "threads": Int % Range(1, None),
        "verbose": Bool,
        "log": Bool,
        "file_buffer_size": Int % Range(1, None),
        "ignore_warnings": Bool,
        "no_parse_seqids": Bool
    },
    parameter_descriptions={
        "threads": "Number of CPU threads. By default, the program will "
                   "auto-detect and use all available virtual cores on the "
                   "machine.",
        "verbose": "Enable more verbose terminal output.",
        "log": "Enable even more verbose terminal output, which is also "
               "written to a file named diamond.log is the current working "
               "directory.",
        "file_buffer_size": "file buffer size in bytes (default=67108864)",
        "ignore_warnings": "Ignore warnings",
        "no_parse_seqids": "Print raw seqids without parsing"
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
        "eggnog_db": "Artifact containing the eggNOG annotation "
                     "database (eggnog.db, eggnog.taxa.db, "
                     "eggnog.taxa.db.traverse.pkl)"
    },
    name="Fetch the databases necessary to run to run the "
         "eggnog-annotate action.",
    description="Downloads eggnog reference database  "
                "using the `download_eggnog_data.py` script from eggNOG. "
                "Here, this script downloads 3 files "
                "(eggnog.db, eggnog.taxa.db, and eggnog.taxa.db.traverse.pkl) "
                "and creates and artifact with them. At least 80 Gb of "
                "storage space is required to run this action. "
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog.eggnog_diamond_search,
    inputs={
        'sequences': SampleData[Contigs] | FeatureData[MAG],
        'diamond_db': ReferenceDB[Diamond],
    },
    parameters={
        'num_cpus': Int,
        'db_in_memory': Bool,
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
    },
    outputs=[
        ('eggnog_hits', SampleData[BLAST6]),
        ('table', FeatureTable[Frequency])
    ],
    name='Run eggNOG search using diamond aligner',
    description="This method performs the steps by which we find our "
                "possible target sequences to annotate using the diamond "
                "search functionality from the eggnog `emapper.py` script",
)

plugin.methods.register_function(
    function=q2_moshpit.eggnog.eggnog_annotate,
    inputs={
        'eggnog_hits': SampleData[BLAST6],
        'eggnog_db': ReferenceDB[Eggnog],
    },
    parameters={
        'db_in_memory': Bool,
    },
    parameter_descriptions={
        'db_in_memory': 'Read eggnog database into memory. The '
                        'eggnog database is very large(>44GB), so this '
                        'option should only be used on clusters or other '
                        'machines with enough memory.',
    },
    outputs=[('ortholog_annotations', FeatureData[NOG])],
    name='Annotate orthologs against eggNOG database',
    description="Apply eggnog mapper to annotate seed orthologs.",
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
    "datasets_version": Str,
    "download": List[Str],
    "download_base_url": Str,
    "download_path": Str,
    "evalue": Float % Range(0, None, inclusive_start=False),
    "force": Bool,
    "limit": Int % Range(1, 20),
    "help": Bool,
    "list_datasets": Bool,
    "long": Bool,
    "metaeuk_parameters": Str,
    "metaeuk_rerun_parameters": Str,
    "miniprot": Bool,
    "offline": Bool,
    "quiet": Bool,
    "restart": Bool,
    "scaffold_composition": Bool,
    "tar": Bool,
    "update_data": Bool,
    "version": Bool,
}
busco_param_descriptions = {
    "mode": "Specify which BUSCO analysis mode to run."
            "Currently only the 'genome' or 'geno' option is supported, "
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
    "datasets_version": "Specify the version of BUSCO datasets, e.g. odb10.",
    "download": "Download dataset. Possible values are a specific dataset "
                "name, 'all', 'prokaryota', 'eukaryota', or 'virus'. "
                "If used together with other command line arguments, "
                "make sure to place this last. Example: '[dataset ...]'.",
    "download_base_url": "Set the url to the remote BUSCO dataset location.",
    "download_path": "Specify local filepath for storing BUSCO dataset "
                     "downloads.",
    "evalue": "E-value cutoff for BLAST searches. "
              "Allowed formats, 0.001 or 1e-03, Default: 1e-03.",
    "force": "Force rewriting of existing files. Must be used when output "
             "files with the provided name already exist.",
    "help": "Show this help message and exit.",
    "limit": "How many candidate regions (contig or transcript) to consider "
             "per BUSCO. Default: 3.",
    "list_datasets": "Print the list of available BUSCO datasets.",
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
    "offline": "To indicate that BUSCO cannot attempt to download files.",
    "quiet": "Disable the info logs, displays only errors.",
    "restart": "Continue a run that had already partially completed.",
    "scaffold_composition": "Writes ACGTN content per scaffold to a file "
                            "`scaffold_composition.txt`.",
    "tar": "Compress some subdirectories with many files to save space.",
    "update_data": "Download and replace with last versions all lineages "
                   "datasets and files necessary to their automated "
                   "selection.",
    "version": "Show this version and exit.",
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
                "information. ",
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
                "prokaryotic genomes.",
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
