# ----------------------------------------------------------------------------
# Copyright (c) 2022-2023, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from qiime2.core.exceptions import ValidationError
import subprocess
import hashlib
from typing import List


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


def colorify(string):
    return "%s%s%s" % ('\033[1;32m', string, "\033[0m")


def compare_md5_hashes(expected_hash: str, path_to_file: str):
    observed_hash = calculate_md5_from_file(path_to_file)
    if observed_hash != expected_hash:
        raise ValidationError(
            "Download error. Data possibly corrupted.\n"
            f"{path_to_file} has an unexpected MD5 hash.\n\n"
            "Expected hash:\n"
            f"{expected_hash}\n\n"
            "Observed hash:\n"
            f"{observed_hash}"
        )


def calculate_md5_from_file(file_path):
    md5_hash = hashlib.md5()
    with open(file_path, 'rb') as f:
        # Read the file in chunks to handle large files
        for chunk in iter(lambda: f.read(4096), b""):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()
