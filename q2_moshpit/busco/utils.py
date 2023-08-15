arguments_with_hyphens = {
        "auto_lineage": "auto-lineage",
        "auto_lineage_euk": "auto-lineage-euk",
        "auto_lineage_prok": "auto-lineage-prok",
        "list_datasets": "list-datasets",
        "update_data": "update-data",
}

def _parse_busco_params(arg_key, arg_val):
    """Creates a list with argument and its value to be consumed by MetaBAT 2.

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

    # If the key is one of 
    if arg_key in arguments_with_hyphens.keys():
        arg_key = arguments_with_hyphens[arg_key]

    if isinstance(arg_val, bool) and arg_val:
        return [f"--{arg_key}"]
    else:
        return [f"--{arg_key}", str(arg_val)]