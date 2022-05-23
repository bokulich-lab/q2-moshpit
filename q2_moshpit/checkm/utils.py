# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from collections import defaultdict
from typing import Mapping, Dict


def _process_checkm_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by CheckM.

    Argument names will be converted to command line parameters by
    appending a '--' prefix and concatenating words separated by a '_',
    e.g.: 'some_parameter_x' -> '--someParameterX'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """
    if isinstance(arg_val, bool) and arg_val:
        return [f'--{arg_key}']
    else:
        return [f'--{arg_key}', str(arg_val)]


def _get_plots_per_sample(
        all_plots: Mapping[str, Mapping[str, str]]
) -> Dict[str, Dict[str, str]]:
    """Converts mapping of different plot types-to-samples-to-plot paths into
        a new mapping of samples-to-plot types-to-plot paths.

    Args:
        all_plots (Mapping[str, Mapping[str, str]]): Dictionary containing plot
            paths per sample per plot type.

    Returns:
        Dict[str, Dict[str, str]]: Dictionary containing a new mapping of plot
            paths per plot type per sample.
    """
    # check if all plot types have the same count of samples
    all_samples = [set(k.keys()) for k in all_plots.values()]
    if not all([x == all_samples[0] for x in all_samples]):
        raise ValueError(
            'All plot types need to have the same set of samples. '
            f'Sample counts were: '
            f'{[len(k.keys()) for k in all_plots.values()]}.'
        )

    plots_per_sample = defaultdict(dict)
    for key, val in all_plots.items():
        for subkey, subval in val.items():
            plots_per_sample[subkey][key] = subval
    return plots_per_sample
