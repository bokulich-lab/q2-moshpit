# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from q2_moshpit._utils import _process_common_input_params
from q2_moshpit.checkm.utils import _process_checkm_arg

from q2_types_genomics.per_sample_data._format import MultiFASTADirectoryFormat


def evaluate_bins(
    output_dir: str, bins: MultiFASTADirectoryFormat,
    reduced_tree: bool = None, unique: int = None, multi: int = None,
    force_domain: bool = None, no_refinement: bool = None,
    individual_markers: bool = None, skip_adj_correction: bool = None,
    skip_pseudogene_correction: bool = None, aai_strain: float = None,
    ignore_thresholds: bool = None, e_value: float = None,
    length: float = None, threads: int = None,
):

    kwargs = {k: v for k, v in locals().items() if k not in ['bins']}
    common_args = _process_common_input_params(
        processing_func=_process_checkm_arg, params=kwargs
    )

    # TODO: remove later
    return common_args

    # return _bin_contigs_metabat(
    #     contigs=contigs, alignment_maps=alignment_maps,
    #     common_args=common_args
    # )
