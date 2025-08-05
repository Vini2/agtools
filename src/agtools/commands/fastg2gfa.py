#!/usr/bin/env python3

import re

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "0.1.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Alpha"


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
            Mapping from node ID to DNA sequence.
        edges : dict
            Mapping from node ID to list of adjacent node IDs.
    """

    segments = {}
    edges = {}

    with open(fastg_file, "r") as f:
        current_node = None
        sequence = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_node and sequence:
                    segments[current_node] = "".join(sequence)
                    sequence = []

                header = line[1:]
                parts = header.split(":")
                node = parts[0].strip("'")
                neighbors = []
                if len(parts) > 1:
                    neighbors = re.split(r"[,\s]+", parts[1])
                current_node = node
                edges[node] = neighbors
            else:
                sequence.append(line)

        if current_node and sequence:
            segments[current_node] = "".join(sequence)

    return segments, edges


def _extract_links(edges: dict, fixed_overlap: int) -> list:
    """
    Convert an edge adjacency dictionary into GFA-style links.

    Parameters
    ----------
    edges : dict
        Dictionary mapping each segment ID to a list of neighboring segment IDs.
        Segment IDs may include orientation suffixes (e.g., `'` to denote reverse).
    fixed_overlap : int
        The overlap length (in base pairs) to include in each GFA link, formatted as "{overlap}M".

    Returns
    -------
    list of tuple
        A list of tuples representing GFA links. Each tuple contains:
        (from_node, from_orientation, to_node, to_orientation, overlap).
    """
    links = []

    for from_node, neighbors in edges.items():
        for raw_to_node in neighbors:
            if not raw_to_node:
                continue

            from_ori = "+" if not from_node.endswith("'") else "-"
            to_ori = "+" if not raw_to_node.endswith("'") else "-"

            from_clean = from_node.strip("'")
            to_clean = raw_to_node.strip("'")

            links.append((from_clean, from_ori, to_clean, to_ori, f"{fixed_overlap}M"))

    return links


def _write_gfa(segments: dict, links: list, output_path: str) -> str:
    """
    Write segments and links to a GFA file.

    Parameters
    ----------
    segments : dict
        Dictionary of segment IDs mapped to their nucleotide sequences.
    links : list of tuple
        List of tuples representing GFA links. Each tuple is in the format:
        (from_node, from_orientation, to_node, to_orientation, overlap).
    output_path : str
        Directory path where the output GFA file should be saved.

    Returns
    -------
    str
        Full path to the written GFA file.
    """

    output_file = f"{output_path}/converted_graph.gfa"
    with open(output_file, "w") as f:
        for seg_id, seq in segments.items():
            f.write(f"S\t{seg_id}\t{seq}\n")
        for from_node, from_ori, to_node, to_ori, overlap in links:
            f.write(f"L\t{from_node}\t{from_ori}\t{to_node}\t{to_ori}\t{overlap}\n")

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
    links = _extract_links(edges, fixed_overlap=k_overlap)
    output_file = _write_gfa(segments, links, gfa_path)

    return output_file
