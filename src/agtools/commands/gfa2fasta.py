#!/usr/bin/env python3

import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def _get_segment_sequences(gfa_file) -> list:
    """
    Extracts segment sequences from a GFA (Graphical Fragment Assembly) file.

    This function reads the GFA file and extracts lines starting with 'S'
    (which represent segments), cleaning the sequence strings to include only
    valid nucleotide characters (G, A, T, C). It returns a list of BioPython
    SeqRecord objects representing each segment.

    Args:
        gfa_file (str): Path to the input GFA file.

    Returns:
        list: A list of SeqRecord objects containing the cleaned segment sequences.
    """

    sequences = []

    with open(gfa_file) as file:
        line = file.readline()

        while line != "":
            if "S" in line:
                strings = line.split("\t")

                record = SeqRecord(
                    Seq(re.sub("[^GATC]", "", str(strings[2]).upper())),
                    id=str(strings[1]),
                    name=str(strings[1]),
                    description="",
                )

                sequences.append(record)

            line = file.readline()

    return sequences


def _write_segment_sequences(sequences, output_path):
    """
    Writes segment sequences to a FASTA file.

    This function saves a list of BioPython SeqRecord objects to a FASTA file
    named 'segments.fasta' in the specified output directory.

    Args:
        sequences (list): A list of SeqRecord objects to write.
        output_path (str): Directory path where the FASTA file will be saved.

    Returns:
        str: Path to the output FASTA file.
    """
    output_file = f"{output_path}/segments.fasta"
    with open(f"{output_file}", "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")


def gfa2fasta(gfa_file, output_path) -> str:
    """
    Converts a GFA file to a FASTA file containing segment sequences.

    This function reads a GFA file, extracts the segment sequences, and writes
    them to a FASTA file in the specified output directory.

    Args:
        gfa_file (str): Path to the input GFA file.
        output_path (str): Directory path where the output FASTA file should be saved.

    Returns:
        str: Path to the generated FASTA file.
    """
    segment_sequences = _get_segment_sequences(gfa_file)
    output_file = _write_segment_sequences(segment_sequences, output_path)

    return output_file
