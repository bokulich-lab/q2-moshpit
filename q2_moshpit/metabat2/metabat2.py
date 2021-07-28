# ----------------------------------------------------------------------------
# Copyright (c) 2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import glob
import os.path
import shutil
import tempfile

from q2_types_genomics.per_sample_data import ContigSequencesDirFmt, BAMDirFmt
from q2_types_genomics.per_sample_data._format import MultiFASTADirectoryFormat

from q2_moshpit._utils import run_command, _process_common_input_params
from q2_moshpit.metabat2.utils import _process_metabat2_arg


def _get_sample_name_from_path(fp):
    return os.path.basename(fp).split('_')[0]


def _assert_samples(contigs_fps, maps_fps) -> dict:
    # TODO: simplify this with a set
    if len(contigs_fps) != len(maps_fps):
        raise Exception('The number of samples must be equal for both, '
                        'contigs and alignment maps. You provided contigs '
                        f'for {len(contigs_fps)} samples and maps for '
                        f'{len(maps_fps)} samples. Please check your inputs '
                        'and try again.')
    contigs_fps, maps_fps = sorted(contigs_fps), sorted(maps_fps)
    contig_samps = [_get_sample_name_from_path(x) for x in contigs_fps]
    map_samps = [_get_sample_name_from_path(x) for x in maps_fps]
    if all([x == y for x, y in zip(contig_samps, map_samps)]):
        return {s: {'contigs': contigs_fps[i],
                    'map': maps_fps[i]}
                for i, s in enumerate(contig_samps)}
    else:
        raise Exception('Contigs and alignment maps should belong to the '
                        'same sample set. You provided contigs for '
                        f'samples: {",".join(contig_samps)} but maps for '
                        f'samples: {",".join(map_samps)}. Please check '
                        'your inputs and try again.')


def _rename_bin(bin_fp, new_location):
    base, _ = os.path.splitext(os.path.basename(bin_fp))
    return os.path.join(new_location, f'{base.replace(".", "")}.fasta')


def bin_contigs_metabat(
    contigs: ContigSequencesDirFmt, alignment_maps: BAMDirFmt,
    min_contig: int = None, max_p: int = None, min_s: int = None,
    max_edges: int = None, p_tnf: int = None, no_add: bool = None,
    min_cv: int = None, min_cv_sum: int = None, min_cls_size: int = None,
    num_threads: int = None, seed: int = None, debug: bool = None,
    verbose: bool = None
) -> MultiFASTADirectoryFormat:
    kwargs = {k: v for k, v in locals().items()
              if k not in ['contigs', 'alignment_maps']}
    common_args = _process_common_input_params(
        processing_func=_process_metabat2_arg, params=kwargs
    )

    contigs_fps = sorted(glob.glob(os.path.join(str(contigs), '*.fa')))
    maps_fps = sorted(glob.glob(os.path.join(str(alignment_maps), '*.bam')))
    sample_set = _assert_samples(contigs_fps, maps_fps)

    result = MultiFASTADirectoryFormat()
    with tempfile.TemporaryDirectory() as tmp:
        for samp, props in sample_set.items():
            # sort alignment map
            sorted_bam = os.path.join(tmp, f'{samp}_alignment_sorted.bam')
            run_command(
                ['samtools', 'sort', props['map'], '-o', sorted_bam],
                verbose=True
            )
            props['map'] = sorted_bam

            # calculate depth
            depth_fp = os.path.join(str(tmp), f'{samp}_depth.txt')
            run_command(
                ['jgi_summarize_bam_contig_depths', '--outputDepth',
                 depth_fp, props['map']],
                verbose=True
            )

            # run metabat2
            bins_dp = os.path.join(tmp, samp)
            bins_prefix = os.path.join(bins_dp, 'bin')
            os.makedirs(bins_dp)
            cmd = ['metabat2', '-i', props['contigs'], '-a', depth_fp,
                   '-o', bins_prefix]
            cmd.extend(common_args)
            run_command(cmd, verbose=True)

            all_bins = glob.glob(os.path.join(bins_dp, '*.fa'))
            new_location = os.path.join(str(result), samp)
            os.makedirs(new_location)
            all_bins_new = [_rename_bin(x, new_location) for x in all_bins]

            for old, new in zip(all_bins, all_bins_new):
                shutil.move(old, new)

    return result
