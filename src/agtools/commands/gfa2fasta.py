#!/usr/bin/env python3

import re

from rsbio_seq import Sequence

from agtools import __version__
from agtools.commands._format_checks import validate_gfa_input
from agtools.commands._output import open_output_file
from agtools.core._io import open_text

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = __version__
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"

LINE_WIDTH = 60


def _get_segment_sequences(gfa_file: str) -> list:
    """
    Extract segment sequences from a GFA (Graphical Fragment Assembly) file.

    This function reads the GFA file and extracts lines starting with 'S'
    (which represent segments), cleaning the sequence strings to include only
    valid nucleotide characters (G, A, T, C). It returns a list of
    ``rsbio_seq.Sequence`` objects representing each segment. The GFA file may
    be plain text or gzip-compressed.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.

    Returns
    -------
    list of rsbio_seq.Sequence
        A list of Sequence objects containing the cleaned segment sequences.
    """

    sequences = []

    with open_text(gfa_file, "rt", buffering=1024 * 1024) as file:
        for line in file:
            if line.startswith("S"):
                strings = line.split("\t")

                sequences.append(
                    Sequence(
                        str(strings[1]),
                        "",
                        re.sub("[^GATC]", "", str(strings[2].strip()).upper()),
                    )
                )

    return sequences


def _write_segment_sequences(sequences: list, output_path: str) -> str:
    """
    Write segment sequences to a FASTA file.

    This function saves a list of ``rsbio_seq.Sequence`` objects to a FASTA
    file at the specified output path. Sequences are wrapped at
    ``LINE_WIDTH`` characters.

    Parameters
    ----------
    sequences : list of rsbio_seq.Sequence
        A list of Sequence objects to write to the FASTA file.
    output_path : str
        Path where the FASTA file will be saved, or `-` for stdout.

    Returns
    -------
    str
        Path to the output FASTA file.
    """
    with open_output_file(output_path) as (output_file, output_handle):
        for record in sequences:
            output_handle.write(f">{record.id}\n")
            sequence = record.seq
            for start in range(0, len(sequence), LINE_WIDTH):
                output_handle.write(f"{sequence[start : start + LINE_WIDTH]}\n")

    return output_file


def gfa2fasta(gfa_file: str, output_path: str) -> str:
    """
    Convert a GFA file to a FASTA file containing segment sequences.

    This function reads a GFA file, extracts the segment sequences, and writes
    them to a FASTA file at the specified output path.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file. May be plain text or gzip-compressed.
    output_path : str
        Path where the output FASTA file should be saved.

    Returns
    -------
    str
        Path to the generated FASTA file.
    """
    validate_gfa_input(gfa_file, "gfa2fasta")
    segment_sequences = _get_segment_sequences(gfa_file)
    output_file = _write_segment_sequences(segment_sequences, output_path)

    return output_file
