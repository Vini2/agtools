#!/usr/bin/env python3

import re

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def _parse_fastg(fastg_file: str) -> tuple:
    """
    Parse a FASTG file and extract segment sequences and edges.

    Parameters
    ----------
    fastg_file : str
        Path to the FASTG file.

    Returns
    -------
    tuple
        segments : dict
            Mapping from segment ID to DNA sequence.
        edges : dict
            Mapping from segment ID to list of adjacent segment IDs.
    """

    segments = {}  # segment_ID -> sequence
    edges = []  # (from_segment, from_orient, to_segment, to_orient)

    current_segment = None
    segment_key = None
    seq_lines = []

    with open(fastg_file, "r") as f:
        for line in f.readlines():
            line = line.strip()

            if line.startswith(">"):

                # Save previous segment
                if segment_key and "'" not in current_segment:
                    segments[segment_key] = "".join(seq_lines)
                seq_lines = []

                # If there are links
                if ":" in line:
                    parts = line[1:].strip().split(":")
                    from_segment = parts[0]

                    current_segment = from_segment
                    # Remove trailing apostrophe if present
                    segment_key = (
                        current_segment[:-1]
                        if "'" in current_segment
                        else current_segment
                    )

                    from_orientation = "+" if "'" not in from_segment else "-"
                    from_segment_for_link = (
                        from_segment[:-1] if "'" in from_segment else from_segment
                    )

                    to_segments = parts[1].split(",")

                    for to_segment in to_segments:
                        to_orientation = "-" if "'" in to_segment else "+"
                        # Remove trailing semicolon if present
                        to_segment_for_link = (
                            to_segment[:-1] if ";" in to_segment else to_segment
                        )
                        # Remove trailing apostrophe if present
                        to_segment_for_link = (
                            to_segment_for_link[:-1]
                            if "'" in to_segment_for_link
                            else to_segment_for_link
                        )
                        edges.append(
                            (
                                from_segment_for_link,
                                from_orientation,
                                to_segment_for_link,
                                to_orientation,
                            )
                        )
                # If no links, just a segment
                else:
                    # Remove trailing semicolon
                    current_segment = line[1:].strip()[:-1]
                    # Remove trailing apostrophe if present
                    segment_key = (
                        current_segment[:-1]
                        if "'" in current_segment
                        else current_segment
                    )

            else:
                seq_lines.append(line)

    return segments, edges


def _write_gfa(
    segments: dict, edges: list, output_path: str, fixed_overlap: int
) -> str:
    """
    Write segments and links to a GFA file.

    Parameters
    ----------
    segments : dict
        Dictionary of segment IDs mapped to their nucleotide sequences.
    edges : list of tuple
        List of tuples representing GFA links. Each tuple is in the format:
        (from_segment, from_orientation, to_segment, to_orientation, overlap).
    output_path : str
        Directory path where the output GFA file should be saved.

    Returns
    -------
    str
        Full path to the written GFA file.
    """

    output_file = f"{output_path}/converted_graph.gfa"
    with open(output_file, "w") as f:
        f.write("H\tVN:Z:1.0\n")
        for segment_id, sequence in segments.items():
            f.write(f"S\t{segment_id}\t{sequence}\n")
        for from_seg, from_orient, to_seg, to_orient in edges:
            f.write(
                f"L\t{from_seg}\t{from_orient}\t{to_seg}\t{to_orient}\t{fixed_overlap}M\n"
            )

    return output_file


def fastg2gfa(fastg_path: str, k_overlap: int, gfa_path: str) -> str:
    """
    Convert a FASTG file to a GFA file format with fixed k-mer overlap.

    Parameters
    ----------
    fastg_path : str
        Path to the input FASTG file.
    k_overlap : int
        Fixed k-mer overlap length to apply to all links (e.g., 41).
    gfa_path : str
        Directory path where the output GFA file will be saved.

    Returns
    -------
    str
        Full path to the generated GFA file.
    """

    segments, edges = _parse_fastg(fastg_path)
    output_file = _write_gfa(segments, edges, gfa_path, fixed_overlap=k_overlap)

    return output_file
