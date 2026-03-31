#!/usr/bin/env python3

from collections import defaultdict

from Bio.Seq import Seq

from agtools import __version__
from agtools.commands._format_checks import validate_gfa_input
from agtools.commands._output import prepare_output_file

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"

__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def reverse_orientation(orient: str) -> str:
    """
    Reverse the orientation symbol used in GFA links.

    Parameters
    ----------
    orient : str
        Orientation symbol, either ``+`` or ``-``.

    Returns
    -------
    str
        Opposite orientation symbol.
    """

    return "+" if orient == "-" else "-"


def reverse_complement(sequence: str) -> str:
    """
    Obtain the reverse complement of a DNA sequence.

    Parameters
    ----------
    sequence : str
        DNA sequence.

    Returns
    -------
    str
        Reverse-complemented DNA sequence.
    """

    return str(Seq(sequence).reverse_complement())


def _get_graph_sequences(
    gfa_file: str,
) -> tuple[defaultdict[str, set[str]], dict[str, str]]:
    """
    Parse a GFA file and extract segment sequences and oriented links.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.

    Returns
    -------
    tuple
        graph_nodes : defaultdict[str, set[str]]
            Mapping of oriented nodes (for example ``seg+``) to neighboring
            oriented nodes.
        sequences : dict[str, str]
            Mapping of segment IDs to nucleotide sequences.
    """

    sequences = {}
    graph_nodes = defaultdict(set)

    with open(gfa_file, "r") as file:
        for line in file:
            if line.startswith("S\t"):
                parts = line.rstrip().split("\t")
                sequences[parts[1]] = parts[2]

            elif line.startswith("L\t"):
                parts = line.rstrip().split("\t")
                from_seg, from_orient = parts[1], parts[2]
                to_seg, to_orient = parts[3], parts[4]

                from_node = f"{from_seg}{from_orient}"
                to_node = f"{to_seg}{to_orient}"
                graph_nodes[from_node].add(to_node)

                rev_from = f"{to_seg}{reverse_orientation(to_orient)}"
                rev_to = f"{from_seg}{reverse_orientation(from_orient)}"
                graph_nodes[rev_from].add(rev_to)

    return graph_nodes, sequences


def _format_fastg_node(segment_id: str, orient: str) -> str:
    """Format an oriented segment ID using the FASTG apostrophe convention."""

    return f"{segment_id}'" if orient == "-" else segment_id


def _write_fastg(
    graph_nodes: defaultdict[str, set[str]],
    sequences: dict[str, str],
    output_path: str,
) -> str:
    """
    Write the sequence graph to FASTG format.

    Parameters
    ----------
    graph_nodes : defaultdict[str, set[str]]
        Mapping of oriented nodes to their outgoing neighbors.
    sequences : dict[str, str]
        Mapping of segment IDs to nucleotide sequences.
    output_path : str
        Path where the output FASTG file should be saved.

    Returns
    -------
    str
        Full path to the written FASTG file.
    """

    output_file = prepare_output_file(output_path)

    with open(output_file, "w") as file:
        for segment_id, sequence in sequences.items():
            for orient in ("+", "-"):
                node = f"{segment_id}{orient}"
                neighbors = sorted(graph_nodes.get(node, set()))
                header = f">{_format_fastg_node(segment_id, orient)}"

                if neighbors:
                    header += ":" + ",".join(
                        _format_fastg_node(neighbor[:-1], neighbor[-1])
                        for neighbor in neighbors
                    )

                header += ";"

                oriented_sequence = sequence
                if orient == "-" and sequence != "*":
                    oriented_sequence = reverse_complement(sequence)

                file.write(f"{header}\n{oriented_sequence}\n")

    return output_file


def gfa2fastg(gfa_file: str, output_path: str) -> str:
    """
    Convert a GFA file to FASTG format.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.
    output_path : str
        Path where the output FASTG file should be saved.

    Returns
    -------
    str
        Full path to the generated FASTG file.
    """

    validate_gfa_input(gfa_file, "gfa2fastg")
    graph_nodes, sequences = _get_graph_sequences(gfa_file)
    output_file = _write_fastg(graph_nodes, sequences, output_path)

    return output_file
