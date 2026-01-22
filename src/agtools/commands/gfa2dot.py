#!/usr/bin/env python3

from agtools.core.unitig_graph import UnitigGraph


def _write_abyss_dot(graph, output_path):
    """
    Write the graph to a DOT file in ABySS-compatible format.

    Parameters
    ----------
    graph : igraph.Graph
        The unitig graph to export. Vertices should have a 'sequence' attribute.
    output_path : str
        Path to the output directory where the DOT file will be written.

    Returns
    -------
    str
        Full path to the generated DOT file.

    References
    ----------
    ABySS File Formats - DOT
    [https://github.com/bcgsc/abyss/wiki/ABySS-File-Formats#dot](https://github.com/bcgsc/abyss/wiki/ABySS-File-Formats#dot)
    """

    output_file = f"{output_path}/graph.gv"

    with open(output_file, "w") as f:

        f.write(f"digraph g {{\n")

        for segment in graph.graph.vs["name"]:
            f.write(f'"{segment}+" [l={graph.segment_lengths[segment]}]\n')
            f.write(f'"{segment}-" [l={graph.segment_lengths[segment]}]\n')

        for link in graph.link_overlap:
            f.write(f'"{link[0]}" -> "{link[1]}" [d=-{graph.link_overlap[link]}]\n')

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
        Path to the output directory where the DOT file will be written.

    Returns
    -------
    str
        Full path to the generated DOT file.
    """

    output_file = f"{output_path}/graph.dot"
    graph.graph.write_dot(output_file)
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
        Path to the directory where the DOT file will be saved.

    Returns
    -------
    str
        Full path to the generated DOT file.
    """

    ug = UnitigGraph.from_gfa(gfa_file)

    output_file = None

    if abyss:
        output_file = _write_abyss_dot(ug, output_path)
    else:
        output_file = _write_dot(ug, output_path)

    return output_file
