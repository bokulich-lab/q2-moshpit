# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import tempfile

import pandas as pd

from q2_assembly._utils import run_commands_with_pipe
from q2_moshpit._utils import run_command
from q2_types.per_sample_sequences import BAMDirFmt


def rpkm(
        df: pd.DataFrame,
        length_col: str = "length",
        read_counts_col: str = "numreads",
) -> pd.Series:
    df['rpk'] = df[read_counts_col] / (df[length_col] / 10**3)
    reads_per_sample = df.groupby("sample_id")[read_counts_col].sum()
    return df['rpk'] * 10**6 / df["sample_id"].map(reads_per_sample)


def tpm(
        df: pd.DataFrame,
        length_col: str = "length",
        read_counts_col: str = "numreads",
) -> pd.Series:
    df['rpk'] = df[read_counts_col] / df[length_col] / 10**3
    rpk_per_sample = df.groupby("sample_id")['rpk'].sum()
    return df['rpk'] / df["sample_id"].map(rpk_per_sample) * 10**6


def _merge_frames(
        coverage_df: pd.DataFrame, lengths_df: pd.DataFrame
) -> pd.DataFrame:
    coverage_summed = coverage_df.groupby(
        ["sample_id", "mag_id"]
    ).sum().reset_index(drop=False)
    coverage_summed = coverage_summed.merge(
        lengths_df, left_on="mag_id", right_index=True
    )
    return coverage_summed


def _calculate_coverage(
        sample_fp: str, sample_id: str, temp_dir: str,
        min_mapq: int, min_query_len: int, min_base_quality: int,
        min_read_len: int, threads: int
) -> pd.DataFrame:
    """
    Calculate the coverage of a sample.

    This function sorts the sample file using samtools, calculates the
    coverage, and then reads the coverage file into a DataFrame.

    Args:
        sample_fp (str): The file path of the sample.
        sample_id (str): The ID of the sample.
        temp_dir (str): The directory to store temporary files.
        threads (int): The number of threads to use.

    Returns:
        pd.DataFrame: A DataFrame containing the coverage information.
    """
    sample_dir = os.path.join(temp_dir, sample_id)
    output_fp = os.path.join(str(temp_dir), f"{sample_id}.bam")
    coverage_fp = os.path.join(str(temp_dir), f"{sample_id}.coverage.tsv")

    os.makedirs(sample_dir)

    # sort the BAM file
    run_commands_with_pipe(
        cmd1=[
            "samtools", "view",
            "-q", str(min_mapq), "-m", str(min_query_len),
            "--threads", str(threads), sample_fp
        ],
        cmd2=[
            "samtools", "sort",
            "-o", output_fp, "--threads", str(threads), sample_fp
        ],
        verbose=True
    )

    # run_command(
    #     ["samtools", "sort", "-o", output_fp,
    #      "--threads", str(threads), sample_fp],
    #     verbose=True
    # )

    # calculate the coverage
    run_command(
        cmd=[
            "samtools", "coverage",
            "-Q", str(min_base_quality), "-l", str(min_read_len),
            "-o", coverage_fp, output_fp
        ],
        verbose=True
    )

    df = pd.read_csv(coverage_fp, sep="\t", index_col=0)
    df["sample_id"] = sample_id
    df["mag_id"] = df.index.map(lambda x: x.split("_", maxsplit=1)[0])

    return df


def estimate_mag_abundance(
        maps: BAMDirFmt,
        mag_lengths: pd.DataFrame,
        metric: str = "rpkm",
        min_mapq: int = 0,
        min_query_len: int = 0,
        min_base_quality: int = 0,
        min_read_len: int = 0,
        threads: int = 1,
) -> pd.DataFrame:
    metric_func = {"rpkm": rpkm, "tpm": tpm}[metric]

    sample_ids_bam = {
        os.path.basename(x).split("_alignment")[0]: x for x
        in glob.glob(os.path.join(str(maps), "*.bam"))
    }

    # calculate coverage for each sample
    with tempfile.TemporaryDirectory() as temp_dir:
        dfs = []
        for sample_id, sample_fp in sample_ids_bam.items():
            dfs.append(
                _calculate_coverage(
                    sample_fp, sample_id, temp_dir, min_mapq, min_query_len,
                    min_base_quality, min_read_len, threads
                )
            )

    coverage_df = pd.concat(dfs)
    coverage_summed = _merge_frames(coverage_df, mag_lengths)
    coverage_summed["abundance"] = metric_func(coverage_summed)

    # transform into a feature table
    feature_table = coverage_summed.pivot(
        index='sample_id', columns='mag_id', values='abundance'
    )
    feature_table.index.name = "sample-id"

    return feature_table
