#!/usr/bin/env python3

from collections import defaultdict
from Bio.Seq import Seq

def reverse_orientation(orient) -> str:
    """
    Reverses the orientation symbol used in GFA links.

    Args:
        orient (str): Orientation symbol, either '+' or '-'.

    Returns:
        str: The opposite orientation symbol.
    """
    return "+" if orient == "-" else "-"

def reverse_complement(sequence) -> str:
    """
    Obtains the reverse complement of a DNA sequence.

    Args:
        seq (str): DNA sequence.

    Returns:
        str: Reverse complement of the input sequence.
    """
    return str(Seq(sequence).reverse_complement())


def _get_graph_sequences(gfa_file) -> tuple[defaultdict, dict, dict, int]:
    """
    Parses a GFA file to extract sequence and graph structure information.

    This function builds a directed graph from the GFA's 'L' (link) lines,
    stores sequences from 'S' (segment) lines, and computes overlaps between nodes.

    Args:
        gfa_file (str): Path to the GFA file.

    Returns:
        tuple: A tuple of:
            - graph_nodes (dict): Mapping of each oriented node to its neighbors.
            - sequences (dict): Mapping of segment IDs to their sequences.
            - overlaps (int): Overlap length of each link.
            - overlap_value (int): Overlap value
    """
    sequences = {}  # segment_id → sequence
    graph_nodes = defaultdict(set)  # oriented node → set of oriented neighbors
    overlaps = {}  # (from_node, to_node) → int
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

def _write_to_fastg(graph_nodes, sequences, output_path) -> str:
    """
    Writes the sequence graph to a FASTG file format.

    Each node is written with its sequence and connections to neighboring nodes.

    Args:
        graph_nodes (dict): Mapping of each oriented node to its neighbors.
        sequences (dict): Mapping of segment IDs to their sequences.
        output_path (str): Directory to write the FASTG file.

    Returns:
        str: Path to the generated FASTG file.
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

def gfa2fastg(gfa_file, output_path) -> tuple[str, int]:
    """
    Converts a GFA file to a FASTG file representing the sequence graph.

    Parses the GFA file, extracts sequences and graph structure, and writes
    them into FASTG format.

    Args:
        gfa_file (str): Path to the GFA file.
        output_path (str): Directory where the FASTG file will be saved.

    Returns:
        str: Path to the created FASTG file.
    """
    graph_nodes, sequences, overlaps, overlap_value = _get_graph_sequences(gfa_file)
    output_file = _write_to_fastg(graph_nodes, sequences, output_path)

    return output_file, overlap_value