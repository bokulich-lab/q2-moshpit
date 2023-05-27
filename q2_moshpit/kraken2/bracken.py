# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os
import re
import subprocess
import tempfile

import pandas as pd

from q2_moshpit._utils import run_command
from q2_types_genomics.kraken2 import (
    Kraken2ReportDirectoryFormat,
    BrackenDBDirectoryFormat,
)


def _run_bracken_one_sample(
        bracken_db: str, kraken2_report_dir: str, tmp_dir: str,
        threshold: int, read_len: int, level: str
) -> pd.DataFrame:
    sample_id = os.path.basename(kraken2_report_dir)
    kraken2_report_fp = os.path.join(
        kraken2_report_dir, f"{sample_id}.report.txt"
    )
    bracken_output_fp = os.path.join(
        tmp_dir, f"{sample_id}.bracken.output.txt"
    )
    bracken_report_fp = os.path.join(
        tmp_dir, f"{sample_id}.bracken.report.txt"
    )
    cmd = [
        "bracken",
        "-d", bracken_db,
        "-i", kraken2_report_fp,
        "-o", bracken_output_fp,
        "-w", bracken_report_fp,
        "-t", str(threshold),
        "-r", str(read_len),
        "-l", level,
    ]
    try:
        run_command(cmd=cmd, verbose=True)
    except subprocess.CalledProcessError as e:
        # TODO: what should be the behaviour when no reads could be classified?
        raise Exception(
            "An error was encountered while running Bracken, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
    bracken_table = pd.read_csv(bracken_output_fp, sep="\t", index_col=0)
    bracken_table["sample_id"] = sample_id

    return bracken_table


def _estimate_bracken(
        kraken_reports: Kraken2ReportDirectoryFormat,
        bracken_db: BrackenDBDirectoryFormat,
        threshold: int,
        read_len: int,
        level: str
) -> pd.DataFrame:
    bracken_tables = []
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            for report_dir in kraken_reports.path.iterdir():
                bracken_table = _run_bracken_one_sample(
                    bracken_db=str(bracken_db),
                    kraken2_report_dir=report_dir,
                    tmp_dir=tmpdir, threshold=threshold,
                    read_len=read_len, level=level
                )
                bracken_tables.append(bracken_table)
        except subprocess.CalledProcessError as e:
            raise Exception(
                "An error was encountered while running Bracken, "
                f"(return code {e.returncode}), please inspect "
                "stdout and stderr to learn more."
            )

    bracken_table = pd.concat(bracken_tables).reset_index()
    bracken_table = bracken_table.pivot(
        index="sample_id", columns="name", values="new_est_reads"
    )

    return bracken_table


def _assert_read_lens_available(
        bracken_db: BrackenDBDirectoryFormat, read_len: int
):
    pattern = r'.+database(\d{2,})mers\.kmer_distrib$'
    lengths = []
    for db in bracken_db.path.iterdir():
        lengths.extend(re.findall(pattern, str(db)))
    lengths = sorted([int(x) for x in lengths])
    if read_len not in lengths:
        raise ValueError(
            f"Provided read length ({read_len}) is not available in the "
            f"Bracken DB. The available values are: "
            f"{', '.join([str(x) for x in lengths])}."
        )


def classify_kraken_bracken(
    ctx, seqs, kraken2_db, bracken_db, threads=1, confidence=0.0,
    minimum_base_quality=0, memory_mapping=False, minimum_hit_groups=2,
    quick=False, report_minimizer_data=False, threshold=0, read_len=100,
    level='S'
):
    bracken_db = bracken_db.view(BrackenDBDirectoryFormat)
    _assert_read_lens_available(bracken_db, read_len)

    kwargs = {
        k: v for k, v in locals().items()
        if k not in ["ctx", "bracken_db", "threshold", "read_len", "level"]
    }

    classify_kraken = ctx.get_action('moshpit', 'classify_kraken')
    k2_reports, k2_outputs = classify_kraken(**kwargs)

    bracken_table = _estimate_bracken(
        kraken_reports=k2_reports.view(Kraken2ReportDirectoryFormat),
        bracken_db=bracken_db, threshold=threshold, read_len=read_len,
        level=level
    )

    # TODO: should this return the original Kraken2 reports?
    bracken_table = ctx.make_artifact(
        "FeatureTable[Frequency]", bracken_table
    )
    return k2_reports, k2_outputs, bracken_table
