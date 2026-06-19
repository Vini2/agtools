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


def gfa2adj(gfa_file: str, delimiter: str, output_path: str) -> str:
    """
    Convert a GFA file into an adjacency matrix and save it to a file.

    This function parses a GFA file to build an undirected graph representation of the
    assembly, computes its adjacency matrix, and writes the matrix to a tab-separated
    file with segment IDs as row and column headers.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.
    output_path : str
        Path where the output adjacency matrix file will be written.

    Returns
    -------
    str
        Path to the generated file containing the adjacency matrix.
    """

    validate_gfa_input(gfa_file, "gfa2adj")
    ug = UnitigGraph.from_gfa(gfa_file)

    adj_df = ug.get_adjacency_matrix(type="pandas")

    separator = "," if delimiter == "comma" else "\t"
    with open_output_file(output_path) as (output_file, output_handle):
        adj_df.to_csv(output_handle, sep=separator)

    return output_file
