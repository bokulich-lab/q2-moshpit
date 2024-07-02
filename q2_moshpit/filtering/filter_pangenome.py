# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os
import shutil
import subprocess
import tarfile
import tempfile
import urllib
import urllib.error
import urllib.request

import skbio
from tqdm import tqdm

from q2_moshpit._utils import run_command

CHUNK_SIZE = 8192
EBI_SERVER_URL = "ftp://ftp.sra.ebi.ac.uk/vol1/analysis/ERZ127/ERZ12792464/hprc-v1.0-pggb.gfa.gz"
ERR_MSG = (
    "Unable to connect to the EBI server. Please try again later. "
    "The error was: {}"
)


def _fetch_and_extract_pangenome(uri: str, dest_dir: str):
    """
    Fetches and extracts the human pangenome GFA file.

    Args:
        uri (str): The URI of the genome to fetch.
        dest_dir (str): The directory where the data will be saved.

    """
    latest_db = os.path.basename(uri)
    dest_fp = os.path.join(dest_dir, latest_db)
    try:
        response = urllib.request.urlopen(uri)
        total_size = int(response.info().get("Content-Length").strip())
        if total_size > 0:
            progress_bar = tqdm(
                desc=f'Downloading the "{latest_db}" database',
                total=total_size,
                unit="B",
                unit_scale=True,
                unit_divisor=1024,
            )

        with open(dest_fp, "wb") as file:
            while True:
                chunk = response.read(CHUNK_SIZE)
                if not chunk:
                    break
                file.write(chunk)
                if total_size > 0:
                    progress_bar.update(len(chunk))
            progress_bar.close() if total_size > 0 else False
    except urllib.error.URLError as e:
        raise Exception(ERR_MSG.format(e))

    print("Download finished. Extracting database files...")
    run_command(["gunzip", dest_fp])
    os.remove(dest_fp)


def _extract_fasta_from_gfa(gfa_fp: str, fasta_fp: str):
    cmd = ["gfatools", "gfa2fa", gfa_fp]
    with open(fasta_fp, 'w') as f_out:
        try:
            subprocess.run(cmd, stdout=f_out)
        except Exception as e:
            raise Exception(
                f"Failed to extract the fasta file from the GFA. "
                f"The error was: {e}"
            )


def _fetch_and_extract_grch38(get_ncbi_genomes: callable, dest_dir: str):
    results = get_ncbi_genomes(
        taxon='Homo sapiens',
        only_reference=True,
        assembly_levels=['chromosome'],
        assembly_source='refseq',
        only_genomic=False
    )
    results.genome_assemblies.export_data(dest_dir)
    shutil.move(
        os.path.join(dest_dir, "dna-sequences.fasta"),
        os.path.join(dest_dir, "grch38.fasta")
    )


def _combine_fasta_files(*fasta_in_fp, fasta_out_fp):
    with open(fasta_out_fp, 'w') as f_out:
        for f_in in fasta_in_fp:
            for seq in skbio.io.read(f_in, format="fasta"):
                seq.sequence = seq.sequence.upper()
                skbio.io.write(seq, format="fasta", into=f_out)


def filter_reads_pangenome(
        ctx, reads, index=None, n_threads=1, sensitivity='sensitive'
):
    get_ncbi_genomes = ctx.get_action("rescript", "get_ncbi_genomes")
    build_index = ctx.get_action("quality_control", "bowtie2_build")
    filter_reads = ctx.get_action("quality_control", "filter_reads")

    with tempfile.TemporaryDirectory() as tmp:
        if index is None:
            print("Reference index was not provided - it will be generated.")
            print("Fetching human pangenome GFA file...")
            _fetch_and_extract_pangenome(EBI_SERVER_URL, tmp)

            print("Fetching human GRCh38 reference genome...")
            _fetch_and_extract_grch38(get_ncbi_genomes, tmp)

            print("Converting pangenome GFA to FASTA...")
            gfa_fp = glob.glob(os.path.join(tmp, "*.gfa"))[0]
            pan_fasta_fp = os.path.join(tmp, "pangenome.fasta")
            _extract_fasta_from_gfa(gfa_fp, pan_fasta_fp)

            print("Generating an index of the combined reference...")
            combined_fasta_fp = os.path.join(tmp, "combined.fasta")
            _combine_fasta_files(
                pan_fasta_fp, os.path.join(tmp, "grch38.fasta"),
                fasta_out_fp=combined_fasta_fp
            )
            combined_reference = ctx.make_artifact(
                "FeatureData[Sequence]", combined_fasta_fp
            )
            index, = build_index(
                sequences=combined_reference, n_threads=n_threads
            )

        print("Filtering reads against the index...")
        filtered_reads, = filter_reads(
            demultiplexed_sequences=reads,
            database=index,
            n_threads=n_threads,
            exclude_seqs=True,
            sensitivity=sensitivity
        )

    return filtered_reads, index


