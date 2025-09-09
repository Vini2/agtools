#!/usr/bin/env python3

from agtools.core.unitig_graph import UnitigGraph
from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def gfa2adj(gfa_file: str, delimiter: str, output_path: str) -> str:
    """
    Convert a GFA file into an adjacency matrix and save it as a TSV file.

    This function parses a GFA file to build an undirected graph representation of the
    assembly, computes its adjacency matrix, and writes the matrix to a tab-separated
    file with segment IDs as row and column headers.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.
    output_path : str
        Directory path where the output adjacency matrix file will be written.

    Returns
    -------
    str
        Path to the generated TSV file containing the adjacency matrix.
    """

    ug = UnitigGraph.from_gfa(gfa_file)

    adj_df = ug.get_adjacency_matrix(type="pandas")

    separator = "," if delimiter == "comma" else "\t"
    output_file = (
        f"{output_path}/adjacency_matrix.csv"
        if delimiter == "comma"
        else f"{output_path}/adjacency_matrix.tsv"
    )
    adj_df.to_csv(output_file, sep=separator)

    return output_file
