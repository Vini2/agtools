#!/usr/bin/env python3

from agtools import __version__
from agtools.commands._format_checks import validate_gfa_input
from agtools.commands._output import open_output_file
from agtools.core.unitig_graph import UnitigGraph

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = __version__
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def _write_abyss_dot(graph, output_path):
    """
    Write the graph to a DOT file in ABySS-compatible format.

    Parameters
    ----------
    graph : igraph.Graph
        The unitig graph to export. Vertices should have a 'sequence' attribute.
    output_path : str
        Path to the output file where the DOT content will be written.

    Returns
    -------
    str
        Full path to the generated DOT file.

    References
    ----------
    ABySS File Formats - DOT
    [https://github.com/bcgsc/abyss/wiki/ABySS-File-Formats#dot](https://github.com/bcgsc/abyss/wiki/ABySS-File-Formats#dot)
    """

    with open_output_file(output_path) as (output_file, f):

        f.write(f"digraph g {{\n")

        for segment in graph.graph.vs["name"]:
            f.write(f'"{segment}+" [l={graph.segment_lengths[segment]}]\n')
            f.write(f'"{segment}-" [l={graph.segment_lengths[segment]}]\n')

        for (
            from_id,
            from_orient,
            to_id,
            to_orient,
        ), overlap in graph.link_overlap.items():
            from_name = graph.segment_names[from_id]
            to_name = graph.segment_names[to_id]
            f.write(
                f'"{from_name}{from_orient}" -> "{to_name}{to_orient}" [d=-{overlap}]\n'
            )

        f.write(f"}}")

    return output_file


def _write_dot(graph, output_path):
    """
    Write the graph to a standard DOT file using igraph's built-in method.

    Parameters
    ----------
    graph : igraph.Graph
        The graph to export.
    output_path : str
        Path to the output file where the DOT content will be written.

    Returns
    -------
    str
        Full path to the generated DOT file.
    """

    with open_output_file(output_path) as (output_file, output_handle):
        graph.graph.write_dot(output_handle)
    return output_file


def gfa2dot(gfa_file, abyss, output_path):
    """
    Convert a GFA file into a DOT graph format.

    This function parses a GFA file into a unitig graph and writes it to a DOT file.
    It supports two DOT formats: standard and ABySS-compatible.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.
    abyss : bool
        If True, output in ABySS-compatible DOT format. Otherwise, use standard DOT.
    output_path : str
        Path to the output file where the DOT file will be saved.

    Returns
    -------
    str
        Full path to the generated DOT file.
    """

    validate_gfa_input(gfa_file, "gfa2dot")
    ug = UnitigGraph.from_gfa(gfa_file)

    output_file = None

    if abyss:
        output_file = _write_abyss_dot(ug, output_path)
    else:
        output_file = _write_dot(ug, output_path)

    return output_file
