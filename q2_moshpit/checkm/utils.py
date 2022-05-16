# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from collections import defaultdict


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


def _get_plots_per_sample(all_plots):
    plots_per_sample = defaultdict(dict)
    for key, val in all_plots.items():
        for subkey, subval in val.items():
            plots_per_sample[subkey][key] = subval
    return plots_per_sample
