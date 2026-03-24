#!/usr/bin/env python3

from pathlib import Path


def prepare_output_file(output_path: str) -> str:
    """
    Ensure parent directories exist for a target output file path.

    Parameters
    ----------
    output_path : str
        Full output file path provided by the user.

    Returns
    -------
    str
        Normalized output file path as a string.
    """
    output_file = Path(output_path).expanduser()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    return str(output_file)
