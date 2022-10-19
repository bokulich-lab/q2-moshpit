# ----------------------------------------------------------------------------
# Copyright (c) 2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from q2_moshpit._utils import _construct_param


def _process_kraken2_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by Kraken 2.

    Argument names will be converted to command line parameters by
    appending a '--' prefix and replacing all '_' with '-',
    e.g.: 'some_parameter_x' -> '--some-parameter-x'.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """
    if isinstance(arg_val, bool) and arg_val:
        return [_construct_param(arg_key)]
    elif not isinstance(arg_val, list):
        return [_construct_param(arg_key), str(arg_val)]
    else:
        raise NotImplementedError(
            f'Parsing arguments of type "{type(arg_val)}" is not supported.'
        )
