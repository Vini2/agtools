#!/usr/bin/env python3

import re

from agtools.core.fasta_parser import FastaParser
from agtools.core.unitig_graph import UnitigGraph


def _write_filtered_graph(
    segments_to_remove: set, parser: FastaParser, gfa_file: str, output_path: str
) -> str:
    """
    Write a cleaned GFA file by excluding lines that involve specified segments.

    This function processes a GFA file line by line, excluding:
    - `S` lines (segments) whose ID is in `segments_to_remove`
    - `L`, `J`, and `C` lines (links, joins, containments) involving any removed segments
    - `P` and `W` lines (paths, walks) that contain any removed segment IDs

    All other lines (including headers or comments) are preserved.

    Parameters
    ----------
    segments_to_remove : set[str]
        Segment IDs that should be removed from the graph.
    parser : FastaParser
        FASTA parser.
    gfa_file : str
        Path to the input GFA file.
    output_path : str
        Directory where the cleaned GFA file will be written.

    Returns
    -------
    str
        Full path to the written cleaned GFA file.
    """

    output_file = f"{output_path}/cleaned_graph.gfa"

    with open(gfa_file, "r") as gfa, open(output_file, "w") as cleaned_gfa:
        for line in gfa:
            if line.startswith("S"):
                parts = line.strip().split("\t")
                seg_id = parts[1]

                if seg_id not in segments_to_remove:
                    if parts[2] == "":
                        parts[2] = str(parser.get_sequence(seg_id))
                        line = "\t".join(parts) + "\n"
                    cleaned_gfa.write(line)

            elif line.startswith("L") or line.startswith("J"):
                parts = line.strip().split("\t")
                from_seg, to_seg = parts[1], parts[3]
                if (
                    from_seg not in segments_to_remove
                    and to_seg not in segments_to_remove
                ):
                    cleaned_gfa.write(line)
            elif line.startswith("C"):
                parts = line.strip().split("\t")
                container_seg, contained_seg = parts[1], parts[3]
                if (
                    container_seg not in segments_to_remove
                    and contained_seg not in segments_to_remove
                ):
                    cleaned_gfa.write(line)
            elif line.startswith("P"):
                parts = line.strip().split("\t")
                seg_ids = parts[2].split(",")
                if all(seg_id not in segments_to_remove for seg_id in seg_ids):
                    cleaned_gfa.write(line)
            elif line.startswith("W"):
                parts = line.strip().split("\t")
                seg_ids = re.split(r"[><]", parts[-1])
                if all(seg_id not in segments_to_remove for seg_id in seg_ids):
                    cleaned_gfa.write(line)
            else:
                cleaned_gfa.write(line)

    return output_file


def clean(gfa_file: str, fasta: str, assembler: str, output_path: str) -> str:
    """
    Clean a GFA file based on segments in a FASTA file.

    This function adds the FASTA sequence to the GFA file if
    missing, removes segments if not present in the FASTA file
    and removes any links, paths, walks, junctions or
    containments containing missing segments.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.
    fasta : int
        Path to the FASTA file.
    assembler : str
        Assembler used to get the assembly
    output_path : str
        Directory where the filtered GFA file will be saved.

    Returns
    -------
    str
        Full path to the cleaned GFA file.
    """

    ug = UnitigGraph.from_gfa(gfa_file)

    # Get parser for fasta file
    parser = FastaParser(fasta, assembler=assembler)

    segments_to_remove = set()

    for segment in ug.segment_names:
        if segment not in parser.index:
            segments_to_remove.add(segment)

    output_file = _write_filtered_graph(
        segments_to_remove, parser, gfa_file, output_path
    )

    return output_file
