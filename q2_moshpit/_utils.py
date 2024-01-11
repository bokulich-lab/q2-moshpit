# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import subprocess
from typing import List

# CONVERT shell colors to the same curses palette
SHELL_COLORS = {
    "wr": '\033[1;37;41m',  # white on red
    "wo": '\033[1;37;43m',  # white on orange
    "wm": '\033[1;37;45m',  # white on magenta
    "wb": '\033[1;37;46m',  # white on blue
    "bw": '\033[1;37;40m',  # black on white
    "lblue": '\033[1;34m',  # light blue
    "lred": '\033[1;31m',  # light red
    "lgreen": '\033[1;32m',  # light green
    "yellow": '\033[1;33m',  # yellow
    "cyan": '\033[36m',  # cyan
    "blue": '\033[34m',  # blue
    "green": '\033[32m',  # green
    "orange": '\033[33m',  # orange
    "red": '\033[31m',  # red
    "magenta": "\033[35m",  # magenta
    "white": "\033[0m",  # white
    None: "\033[0m",  # end
}


def run_command(cmd, env=None, verbose=True, pipe=False, **kwargs):
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')

    if pipe:
        result = subprocess.run(
            cmd, env=env, check=True, capture_output=True, text=True
        )
        return result

    if env:
        subprocess.run(cmd, env=env, check=True, **kwargs)
    else:
        subprocess.run(cmd, check=True, **kwargs)


def _construct_param(arg_name):
    """Converts argument name into a command line parameter."""
    return f'--{arg_name.replace("_", "-")}'


def _process_common_input_params(processing_func, params: dict) -> List[str]:
    """Converts provided arguments and their values.

    Conversion is entirely dependent on the passed 'processing_func'
    that processes individual arguments. The output is a list of
    parameters with their values that can be directly passed to the
    respective command.

    Arguments without any value are skipped.
    Any other argument is processed using the 'processing_func' and
    appended to the final list.

    Args:
        processing_func: Function to be used for formatting a single argument.
        params (dict): Dictionary of parameter: value pairs to be processed.

    Returns:
        processed_args (list): List of processed arguments and their values.

    """
    processed_args = []
    for arg_key, arg_val in params.items():
        # This if condition excludes arguments which are falsy
        # (False, None, "", []), except for integers and floats.
        if (  # noqa: E721
            type(arg_val) == int or
            type(arg_val) == float or
            arg_val
        ):
            processed_args.extend(processing_func(arg_key, arg_val))
        else:
            continue

    return processed_args


def colorify(string, color):
    return "%s%s%s" % (SHELL_COLORS[color], string, SHELL_COLORS[None])
