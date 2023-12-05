# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import subprocess
import os
import tempfile
import qiime2.util
# import pdb here for debugging purposes
import pandas as pd
from typing import Union
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
from q2_types_genomics.genome_data import SeedOrthologDirFmt, OrthologFileFmt
from .._utils import run_command, _process_common_input_params
from ._utils import _parse_build_diamond_db_params
from q2_types_genomics.reference_db import EggnogRefDirFmt
from q2_types.feature_data import (
    DNAFASTAFormat, ProteinSequencesDirectoryFormat
)
from q2_types_genomics.reference_db import (
    DiamondDatabaseDirFmt, NCBITaxonomyDirFmt
)
from q2_types_genomics.feature_data import (
    OrthologAnnotationDirFmt, MAGSequencesDirFmt
)


def eggnog_diamond_search(
        sequences: Union[ContigSequencesDirFmt, MAGSequencesDirFmt],
        diamond_db: DiamondDatabaseDirFmt,
        num_cpus: int = 1, db_in_memory: bool = False
) -> (SeedOrthologDirFmt, pd.DataFrame):

    diamond_db_fp = os.path.join(str(diamond_db), 'ref_db.dmnd')
    temp = tempfile.TemporaryDirectory()

    # run analysis
    if isinstance(sequences, ContigSequencesDirFmt):
        for relpath, obj_path in sequences.sequences.iter_views(
                DNAFASTAFormat):
            sample_id = str(relpath).rsplit(r'_', 1)[0]
            _diamond_search_runner(
                input_path=obj_path, diamond_db=diamond_db_fp,
                sample_label=sample_id, output_loc=temp.name,
                num_cpus=num_cpus, db_in_memory=db_in_memory
            )
    elif isinstance(sequences, MAGSequencesDirFmt):
        for mag_fp in glob.glob(f'{sequences.path}/*.fa*'):
            sample_id = os.path.splitext(os.path.basename(mag_fp))[0]
            _diamond_search_runner(
                input_path=mag_fp, diamond_db=diamond_db_fp,
                sample_label=sample_id, output_loc=temp.name,
                num_cpus=num_cpus, db_in_memory=db_in_memory
            )

    result = SeedOrthologDirFmt()
    ortholog_fps = [
        os.path.basename(x) for x
        in glob.glob(f'{temp.name}/*.seed_orthologs')
    ]
    for item in ortholog_fps:
        qiime2.util.duplicate(
            os.path.join(temp.name, item), os.path.join(result.path, item)
        )

    ft = _eggnog_feature_table(result)

    return result, ft


def _eggnog_feature_table(seed_orthologs: SeedOrthologDirFmt) -> pd.DataFrame:
    per_sample_counts = []

    for sample_path, obj in seed_orthologs.seed_orthologs.iter_views(
            OrthologFileFmt):
        # TODO: put filename to sample name logic on OrthologFileFmt object
        sample_name = str(sample_path).replace('.emapper.seed_orthologs', '')
        sample_df = obj.view(pd.DataFrame)
        sample_feature_counts = sample_df.value_counts('sseqid')
        sample_feature_counts.name = str(sample_name)
        per_sample_counts.append(sample_feature_counts)
    df = pd.DataFrame(per_sample_counts)
    df.fillna(0, inplace=True)
    df.columns = df.columns.astype('str')

    return df


def _diamond_search_runner(input_path, diamond_db, sample_label, output_loc,
                           num_cpus, db_in_memory):

    cmds = ['emapper.py', '-i', str(input_path), '-o', sample_label,
            '-m', 'diamond', '--no_annot', '--dmnd_db', str(diamond_db),
            '--itype', 'metagenome', '--output_dir', output_loc, '--cpu',
            str(num_cpus)]
    if db_in_memory:
        cmds.append('--dbmem')

    subprocess.run(cmds, check=True)


def eggnog_annotate(eggnog_hits: SeedOrthologDirFmt,
                    eggnog_db: EggnogRefDirFmt,
                    db_in_memory: bool = False) -> OrthologAnnotationDirFmt:

    eggnog_db_fp = eggnog_db.path

    result = OrthologAnnotationDirFmt()

    # run analysis
    for relpath, obj_path in eggnog_hits.seed_orthologs.iter_views(
            OrthologFileFmt):
        sample_label = str(relpath).rsplit(r'.', 2)[0]

        _annotate_seed_orthologs_runner(seed_ortholog=obj_path,
                                        eggnog_db=eggnog_db_fp,
                                        sample_label=sample_label,
                                        output_loc=result,
                                        db_in_memory=db_in_memory)

    return result


def _annotate_seed_orthologs_runner(seed_ortholog, eggnog_db, sample_label,
                                    output_loc, db_in_memory):

    # at this point instead of being able to specify the type of target
    # orthologs, we want to annotate _all_.

    cmds = ['emapper.py', '-m', 'no_search', '--annotate_hits_table',
            str(seed_ortholog), '--data_dir', str(eggnog_db),
            '-o', str(sample_label), '--output_dir', str(output_loc)]
    if db_in_memory:
        cmds.append('--dbmem')

    subprocess.run(cmds, check=True)


def build_diamond_db(
        sequences: ProteinSequencesDirectoryFormat,
        taxonomy_data: NCBITaxonomyDirFmt = None,
        threads: int = None,
        verbose: bool = False,
        log: bool = False,
        file_buffer_size: int = 67108864,
        ignore_warnings: bool = False,
        no_parse_seqids: bool = False
        ) -> DiamondDatabaseDirFmt:
    '''
    Builds diamond database from protein reference database file in FASTA
    format.
    '''
    # Process input parameters
    kwargs = {}
    for key, value in locals().items():
        if key not in ["sequences", "taxonomy_data", "kwargs"]:
            kwargs[key] = value

    # Add paths to taxonomy data if provided
    if taxonomy_data is not None:
        kwargs["taxonmap"] = os.path.join(
            str(taxonomy_data), "prot.accession2taxid.gz"
            )
        kwargs["taxonnodes"] = os.path.join(str(taxonomy_data), "nodes.dmp")
        kwargs["taxonnames"] = os.path.join(str(taxonomy_data), "names.dmp")

    # Filter out all kwargs that are falsy (except 0 and 0.0)
    parsed_args = _process_common_input_params(
        processing_func=_parse_build_diamond_db_params, params=kwargs
    )

    # Instantiate output object
    diamond_db = DiamondDatabaseDirFmt()

    # Define path to in/output file
    path_in = os.path.join(str(sequences), "protein-sequences.fasta")
    path_out = os.path.join(str(diamond_db), "ref_db.dmnd")

    # Run diamond makedb
    cmd = [
        "diamond", "makedb",
        "--in", f"{path_in}",
        "--db", f"{path_out}"
    ]
    cmd.extend(parsed_args)
    run_command(cmd)

    # Return output artifact
    return diamond_db


def fetch_eggnog_db() -> EggnogRefDirFmt:
    """
    Downloads eggnog reference database using the
    `download_eggnog_data.py` script from eggNOG. Here, this
    script downloads 3 files amounting to 47Gb in total.
    """

    # Initialize output objects
    eggnog_db = EggnogRefDirFmt()

    # Define command.
    # Meaning of flags:
    # y: Answer yest to all prompts thrown by download_eggnog_data.py
    # D: Do not download the Diamond database
    # data_dir: location where to save downloads
    cmd = [
        "download_eggnog_data.py", "-y", "-D",
        "--data_dir", str(eggnog_db.path)
    ]
    run_command(cmd)

    # Return objects
    return eggnog_db
