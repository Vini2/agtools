#!/usr/bin/env python3

import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO

STDOUT_PATH = "-"


def prepare_output_file(output_path: str) -> str:
    """
    Ensure parent directories exist for a target output file path.

    Parameters
    ----------
    output_path : str
        Full output file path provided by the user. A single hyphen (`-`)
        routes output to stdout.

    Returns
    -------
    str
        Normalized output file path as a string, or `-` for stdout.
    """
    if output_path == STDOUT_PATH:
        return STDOUT_PATH

    output_file = Path(output_path).expanduser()
    output_file.parent.mkdir(parents=True, exist_ok=True)
    return str(output_file)


@contextmanager
def open_output_file(output_path: str) -> Iterator[tuple[str, TextIO]]:
    """
    Open an output destination for writing text.

    Parameters
    ----------
    output_path : str
        Full output file path provided by the user. A single hyphen (`-`)
        routes output to stdout.

    Yields
    ------
    tuple[str, TextIO]
        The normalized output identifier and an open writable text handle.
    """
    output_file = prepare_output_file(output_path)

    if output_file == STDOUT_PATH:
        yield output_file, sys.stdout
        return

    with open(output_file, "w") as output_handle:
        yield output_file, output_handle
