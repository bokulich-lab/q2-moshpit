# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import subprocess
from copy import deepcopy
from typing import Union, Optional

import pandas as pd
from q2_types.per_sample_sequences import (
    SingleLanePerSamplePairedEndFastqDirFmt,
    SingleLanePerSampleSingleEndFastqDirFmt
)
from q2_types.feature_data import DNAFASTAFormat

from q2_moshpit._utils import run_command, _process_common_input_params
from q2_moshpit.kraken2.utils import _process_kraken2_arg
from q2_types_genomics.feature_data import MAGSequencesDirFmt
from q2_types_genomics.per_sample_data import ContigSequencesDirFmt
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


def _classify_kraken2(
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
        if manifest is not None:
            # we got reads - use the manifest
            iterate_over = manifest.iterrows()
            path_function = get_paths_for_reads
        else:
            # we got contigs or MAGs
            def view_to_paths(arg):
                relpath, _ = arg
                return relpath.stem, str(seqs.path / relpath)

            iterate_over = map(
                view_to_paths, seqs.sequences.iter_views(DNAFASTAFormat)
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


def classify_kraken2(
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
              if k not in ["seqs", "kraken2_db"]}
    common_args = _process_common_input_params(
        processing_func=_process_kraken2_arg, params=kwargs
    )
    common_args.extend(["--db", str(kraken2_db.path)])
    return _classify_kraken2(seqs, common_args)
