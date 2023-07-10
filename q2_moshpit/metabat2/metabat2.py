# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
from uuid import uuid4

import os.path
import shutil
import tempfile
from copy import deepcopy

import skbio.io
from q2_types.feature_data import DNAIterator

from q2_types_genomics.per_sample_data import ContigSequencesDirFmt, BAMDirFmt
from q2_types_genomics.per_sample_data._format import MultiFASTADirectoryFormat

from q2_moshpit._utils import run_command, _process_common_input_params
from q2_moshpit.metabat2.utils import _process_metabat2_arg


def _get_sample_name_from_path(fp, suffix):
    return os.path.basename(fp).rsplit(suffix, maxsplit=1)[0]


def _assert_samples(contigs_fps, maps_fps) -> dict:
    contigs_fps, maps_fps = sorted(contigs_fps), sorted(maps_fps)
    contig_samps = [_get_sample_name_from_path(x, '_contigs.fa')
                    for x in contigs_fps]
    map_samps = [_get_sample_name_from_path(x, '_alignment.bam')
                 for x in maps_fps]
    if set(contig_samps) != set(map_samps):
        raise Exception('Contigs and alignment maps should belong to the '
                        'same sample set. You provided contigs for '
                        f'samples: {",".join(contig_samps)} but maps for '
                        f'samples: {",".join(map_samps)}. Please check '
                        'your inputs and try again.')
    return {s: {'contigs': contigs_fps[i], 'map': maps_fps[i]}
            for i, s in enumerate(contig_samps)}


def _sort_bams(samp_name, samp_props, loc):
    sorted_bam = os.path.join(loc, f'{samp_name}_alignment_sorted.bam')
    new_props = deepcopy(samp_props)
    run_command(
        ['samtools', 'sort', new_props['map'], '-o', sorted_bam],
        verbose=True
    )
    new_props['map'] = sorted_bam
    return new_props


def _estimate_depth(samp_name, samp_props, loc):
    depth_fp = os.path.join(str(loc), f'{samp_name}_depth.txt')
    run_command(
        ['jgi_summarize_bam_contig_depths', '--outputDepth',
         depth_fp, samp_props['map']],
        verbose=True
    )
    return depth_fp


def _run_metabat2(samp_name, samp_props, loc, depth_fp, common_args):
    bins_dp = os.path.join(loc, samp_name)
    bins_prefix = os.path.join(bins_dp, 'bin')
    os.makedirs(bins_dp)
    cmd = ['metabat2', '-i', samp_props['contigs'], '-a', depth_fp,
           '-o', bins_prefix]
    cmd.extend(common_args)
    run_command(cmd, verbose=True)
    return bins_dp


def _rename_bin(new_location):
    uuid = uuid4()
    return os.path.join(new_location, f'{uuid}.fasta')


def _process_sample(samp_name, samp_props, common_args, result_loc):
    with tempfile.TemporaryDirectory() as tmp:
        # sort alignment map
        props = _sort_bams(samp_name, samp_props, tmp)

        # calculate depth
        depth_fp = _estimate_depth(samp_name, props, tmp)

        # run metabat2
        bins_dp = _run_metabat2(
            samp_name, props, tmp, depth_fp, common_args
        )

        # rename using UUID v4
        all_bins = glob.glob(os.path.join(bins_dp, '*.fa'))
        new_location = os.path.join(str(result_loc), samp_name)
        os.makedirs(new_location)
        all_bins_new = [_rename_bin(new_location) for x in all_bins]

        for old, new in zip(all_bins, all_bins_new):
            shutil.move(old, new)


def _generate_contig_map(
        bins: MultiFASTADirectoryFormat
) -> dict:
    contig_map = {}
    for seq in bins.sequences.iter_views(DNAIterator):
        mag_id = os.path.splitext(os.path.basename(seq[0]))[0]
        seqs = skbio.read(
            os.path.join(str(bins), str(seq[0])),
            format='fasta', verify=False
        )
        contigs = [x.metadata['id'] for x in seqs]
        contig_map[mag_id] = contigs
    return contig_map


def _bin_contigs_metabat(
        contigs: ContigSequencesDirFmt, alignment_maps: BAMDirFmt,
        common_args: list
) -> (MultiFASTADirectoryFormat, dict):
    contigs_fps = sorted(glob.glob(os.path.join(str(contigs), '*.fa')))
    maps_fps = sorted(glob.glob(os.path.join(str(alignment_maps), '*.bam')))
    sample_set = _assert_samples(contigs_fps, maps_fps)

    bins = MultiFASTADirectoryFormat()
    for samp, props in sample_set.items():
        _process_sample(samp, props, common_args, str(bins))

    contig_map = _generate_contig_map(bins)

    return bins, contig_map


def bin_contigs_metabat(
    contigs: ContigSequencesDirFmt, alignment_maps: BAMDirFmt,
    min_contig: int = None, max_p: int = None, min_s: int = None,
    max_edges: int = None, p_tnf: int = None, no_add: bool = None,
    min_cv: int = None, min_cv_sum: int = None, min_cls_size: int = None,
    num_threads: int = None, seed: int = None, debug: bool = None,
    verbose: bool = None
) -> (MultiFASTADirectoryFormat, dict):

    kwargs = {k: v for k, v in locals().items()
              if k not in ['contigs', 'alignment_maps']}
    common_args = _process_common_input_params(
        processing_func=_process_metabat2_arg, params=kwargs
    )

    return _bin_contigs_metabat(
        contigs=contigs, alignment_maps=alignment_maps,
        common_args=common_args
    )
