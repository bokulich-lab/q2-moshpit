# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import subprocess
from copy import deepcopy
from typing import Union, Optional

import pandas as pd
from q2_types.per_sample_sequences import (
    SequencesWithQuality,
    PairedEndSequencesWithQuality,
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt
)
from q2_types.sample_data import SampleData
from q2_types.feature_data import FeatureData

from q2_moshpit._utils import run_command, _process_common_input_params
from q2_moshpit.kraken2.utils import _process_kraken2_arg
from q2_types_genomics.feature_data import MAGSequencesDirFmt, MAG
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt, Contigs
from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    Kraken2OutputDirectoryFormat,
    Kraken2DBDirectoryFormat,
)


def _get_seq_paths(df_index, df_row, df_columns):
    if "reverse" in df_columns:
        _sample, fn = df_index, df_row.tolist()
    else:
        _sample, fn = df_index, [df_row["forward"]]
    return _sample, fn


def _construct_output_paths(
        _sample, kraken2_outputs_dir, kraken2_reports_dir
):
    report_fp = os.path.join(
        kraken2_reports_dir.path, f"{_sample}.report.txt"
    )
    output_fp = os.path.join(
        kraken2_outputs_dir.path, f"{_sample}.output.txt"
    )
    return output_fp, report_fp


def classify_kraken2(
        ctx,
        seqs,
        kraken2_db,
        threads=1,
        confidence=0.0,
        minimum_base_quality=0,
        memory_mapping=False,
        minimum_hit_groups=2,
        quick=False,
        report_minimizer_data=False,
        num_partitions=None
):
    kwargs = {k: v for k, v in locals().items()
              if k not in ["seqs", "kraken2_db", "ctx", "num_partitions"]}

    _classify_kraken2 = ctx.get_action("moshpit", "_classify_kraken2")
    collate_kraken2_reports = ctx.get_action("moshpit",
                                             "collate_kraken2_reports")
    collate_kraken2_outputs = ctx.get_action("moshpit",
                                             "collate_kraken2_outputs")

    if seqs.type <= SampleData[SequencesWithQuality]:
        partition_method = ctx.get_action("demux", "partition_samples_single")
    elif seqs.type <= SampleData[PairedEndSequencesWithQuality]:
        partition_method = ctx.get_action("demux", "partition_samples_paired")
    elif seqs.type <= SampleData[Contigs]:
        partition_method = ctx.get_action("assembly", "partition_contigs")
    # FeatureData[MAG] is not parallelized
    elif seqs.type <= FeatureData[MAG]:
        kraken2_reports, kraken2_outputs = \
                _classify_kraken2(seqs, kraken2_db, **kwargs)
        return kraken2_reports, kraken2_outputs
    else:
        raise NotImplementedError()

    (partitioned_seqs,) = partition_method(seqs, num_partitions)

    kraken2_reports = []
    kraken2_outputs = []
    for seq in partitioned_seqs.values():
        (kraken2_report, kraken2_output) = _classify_kraken2(
                seq, kraken2_db, **kwargs)
        kraken2_reports.append(kraken2_report)
        kraken2_outputs.append(kraken2_output)

    (collated_kraken2_reports,) = collate_kraken2_reports(kraken2_reports)
    (collated_kraken2_outputs,) = collate_kraken2_outputs(kraken2_outputs)

    return collated_kraken2_reports, collated_kraken2_outputs


def _classify_kraken2(
        seqs: Union[
            SingleLanePerSamplePairedEndFastqDirFmt,
            SingleLanePerSampleSingleEndFastqDirFmt,
            ContigSequencesDirFmt,
            MAGSequencesDirFmt,
        ],
        kraken2_db: Kraken2DBDirectoryFormat,
        threads: int = 1,
        confidence: float = 0.0,
        minimum_base_quality: int = 0,
        memory_mapping: bool = False,
        minimum_hit_groups: int = 2,
        quick: bool = False,
        report_minimizer_data: bool = False
) -> (
        Kraken2ReportDirectoryFormat,
        Kraken2OutputDirectoryFormat,
):
    kwargs = {k: v for k, v in locals().items()
              if k not in ["seqs", "kraken2_db", "ctx"]}

    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    common_args.extend(["--db", str(kraken2_db.path)])
    return classify_kraken2_helper(seqs, common_args)


def classify_kraken2_helper(
        seqs, common_args
) -> (Kraken2ReportDirectoryFormat, Kraken2OutputDirectoryFormat):
    if isinstance(seqs, (MAGSequencesDirFmt, ContigSequencesDirFmt)):
        manifest = None
    else:
        manifest: Optional[pd.DataFrame] = seqs.manifest.view(pd.DataFrame)

    base_cmd = ["kraken2", *common_args]
    if manifest is not None and "reverse" in manifest.columns:
        base_cmd.append("--paired")

    kraken2_reports_dir = Kraken2ReportDirectoryFormat()
    kraken2_outputs_dir = Kraken2OutputDirectoryFormat()

    def get_paths_for_reads(index, row):
        return _get_seq_paths(index, row, list(manifest.columns))

    def get_paths_for_mags(mag_id, fp):
        return mag_id, [fp]

    def get_paths_for_contigs(contig_id, fp):
        # HACK: remove after adding manifest or other solution, see
        # https://github.com/bokulich-lab/q2-types-genomics/issues/56
        return contig_id.rstrip('_contigs'), [fp]

    try:
        if manifest is not None:  # we got reads - use the manifest
            iterate_over = manifest.iterrows()
            path_function = get_paths_for_reads
        else:
            iterate_over = (
                (os.path.basename(fp).split(".")[0], fp)
                for fp in sorted(glob.glob(os.path.join(seqs.path, "*.fasta")))
            )
            if type(seqs) is MAGSequencesDirFmt:
                path_function = get_paths_for_mags
            elif type(seqs) is ContigSequencesDirFmt:
                path_function = get_paths_for_contigs

        for args in iterate_over:
            _sample, fn = path_function(*args)
            output_fp, report_fp = _construct_output_paths(
                _sample, kraken2_outputs_dir, kraken2_reports_dir
            )
            cmd = deepcopy(base_cmd)
            cmd.extend(
                ["--report", report_fp, "--output", output_fp, *fn]
            )
            run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running Kraken 2, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )

    return kraken2_reports_dir, kraken2_outputs_dir
