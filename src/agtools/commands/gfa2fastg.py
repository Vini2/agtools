#!/usr/bin/env python3

from collections import defaultdict

from Bio.Seq import Seq

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "0.1.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Alpha"


def reverse_orientation(orient: str) -> str:
    """
    Reverse the orientation symbol used in GFA links.

    Parameters
    ----------
    orient : str
        Orientation symbol, either '+' or '-'.

    Returns
    -------
    str
        The opposite orientation symbol ('-' if input is '+', and vice versa).
    """

    return "+" if orient == "-" else "-"


def reverse_complement(sequence: str) -> str:
    """
    Obtain the reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence consisting of A, T, G, and C characters.

    Returns
    -------
    str
        Reverse complement of the input DNA sequence.
    """

    return str(Seq(sequence).reverse_complement())


def _get_graph_sequences(gfa_file: str) -> tuple[defaultdict, dict, dict, int]:
    """
    Parse a GFA file to extract sequence and graph structure information.

    This function reads a GFA file and constructs a directed graph from 'L' (link) lines,
    stores sequences from 'S' (segment) lines, and extracts overlap information for each link.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.

    Returns
    -------
    tuple
        A tuple containing:

        - graph_nodes : dict
            Dictionary mapping each oriented node (e.g., '1+', '2-') to a list of neighboring nodes.

        - sequences : dict
            Dictionary mapping segment IDs to their nucleotide sequences.

        - overlaps : int
            Overlap length extracted from the GFA links.

        - overlap_value : int
            The fixed or computed overlap value.
    """

    sequences = {}  # segment_id -> sequence
    graph_nodes = defaultdict(set)  # oriented node -> set of oriented neighbors
    overlaps = {}  # (from_node, to_node) -> int
    overlap_value = 0

    with open(gfa_file) as f:
        for line in f:
            if line.startswith("S"):
                parts = line.strip().split("\t")
                seg_id, seq = parts[1], parts[2]
                sequences[seg_id] = seq

            elif line.startswith("L"):
                parts = line.strip().split("\t")
                from_seg, from_orient = parts[1], parts[2]
                to_seg, to_orient = parts[3], parts[4]
                overlap_value = int(parts[5][:-1])  # Remove trailing M

                from_node = f"{from_seg}{from_orient}"
                to_node = f"{to_seg}{to_orient}"
                graph_nodes[from_node].add(to_node)
                overlaps[(from_node, to_node)] = overlap_value

                # Add reverse link
                rev_from = f"{to_seg}{reverse_orientation(to_orient)}"
                rev_to = f"{from_seg}{reverse_orientation(from_orient)}"
                graph_nodes[rev_from].add(rev_to)
                overlaps[(rev_from, rev_to)] = overlap_value

    return graph_nodes, sequences, overlaps, overlap_value


def _write_to_fastg(graph_nodes: dict, sequences: dict, output_path: str) -> str:
    """
    Write the sequence graph to a FASTG file format.

    Each node is written along with its sequence and a list of connections
    to its neighboring nodes.

    Parameters
    ----------
    graph_nodes : dict
        Dictionary mapping each oriented node (e.g., '1+', '2-') to its neighboring nodes.

    sequences : dict
        Dictionary mapping segment IDs to their nucleotide sequences.

    output_path : str
        Directory path where the FASTG file will be saved.

    Returns
    -------
    str
        Path to the generated FASTG file.
    """
    output_file = f"{output_path}/converted_graph.fastg"
    with open(output_file, "w") as out:
        written = set()

        for seg_id in sequences:
            for orient in ("+", "-"):
                node = f"{seg_id}{orient}"
                if node in written:
                    continue
                written.add(node)

                seq = sequences[seg_id]
                if orient == "-":
                    seq = reverse_complement(seq)

                header = f">{node}"
                neighbors = sorted(graph_nodes.get(node, []))
                if neighbors:
                    header += ":" + ",".join(neighbors)
                out.write(f"{header}\n{seq}\n")

    return output_file


def gfa2fastg(gfa_file: str, output_path: str) -> tuple[str, int]:
    """
    Convert a GFA file to a FASTG file representing the sequence graph.

    This function parses the GFA file to extract sequences and the graph structure,
    then writes the graph into FASTG format.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.

    output_path : str
        Directory path where the FASTG file will be saved.

    Returns
    -------
    str
        Path to the generated FASTG file.
    """
    graph_nodes, sequences, overlaps, overlap_value = _get_graph_sequences(gfa_file)
    output_file = _write_to_fastg(graph_nodes, sequences, output_path)

    return output_file, overlap_value
