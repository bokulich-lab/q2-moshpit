# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

def _process_kaiju_arg(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by Kaiju.

    Argument names will be converted to command line parameters by
    appending a '-' prefix.

    Args:
        arg_key (str): Argument name.
        arg_val: Argument value.

    Returns:
        [converted_arg, arg_value]: List containing a prepared command line
            parameter and, optionally, its value.
    """
    if arg_key == "x":
        return ["-x"] if arg_val else ["-X"]
    elif arg_key == "evalue":
        return ["-E", str(arg_val)]
    elif not isinstance(arg_val, list):
        return [f"-{arg_key}", str(arg_val)]
    else:
        raise NotImplementedError(
            f'Parsing arguments of type "{type(arg_val)}" is not supported.'
        )
