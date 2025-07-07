#!/usr/bin/env python3

import re 

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def _get_segment_sequences(gfa_file) -> list:

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
    output_file = f"{output_path}/segments.fasta"
    with open(f"{output_file}", "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def gfa2fasta(gfa_file, output_path) -> str:

    segment_sequences = _get_segment_sequences(gfa_file)
    output_file = _write_segment_sequences(segment_sequences, output_path)

    return output_file
