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


def _write_stats_file(gfa_file: str, stats: dict, output_path: str) -> str:
    """
    Write the statistics to a file.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.
    stats : dict
        Dictionary containing various computed graph statistics.
    output_path : str
        Directory path where the output statistics file will be saved.

    Returns
    -------
    str
        Path to the written statistics file.
    """
    output_file = f"{output_path}/graph_stats.txt"

    with open(output_file, "w") as f:
        # Write basic graph statistics
        f.write(f"Basic graph statistics for {gfa_file}:\n")
        f.write(f"Number of segments: {stats['nsegments']}\n")
        f.write(f"Number of links: {stats['nlinks']}\n")
        f.write(f"Number of self-loops: {stats['nloops']}\n")
        f.write(f"Number of connected components: {stats['ncomponents']}\n")
        f.write(f"Average node degree: {stats['average_node_degree']}\n")
        f.write("\n")
        # Write sequence-based statistics
        f.write(f"Sequence-based statistics for {gfa_file}:\n")
        f.write(f"Total length of segments: {stats['total_length']} bp\n")
        f.write(f"Average segment length: {stats['average_segment_length']} bp\n")
        f.write(f"N50: {stats['n50']} bp\n")
        f.write(f"L50: {stats['l50']} segment(s)\n")
        f.write(f"GC content: {stats['gc_content']:.2%}")

    return output_file


def stats(gfa_file: str, output_path: str) -> str:
    """
    Compute and write summary statistics for an assembly graph in GFA format.

    This function parses the given GFA file using a UnitigGraph object,
    calculates a variety of assembly and graph-level statistics, and writes
    the results to a file in the specified output directory.

    Parameters:
    ----------
    gfa_file : str
        Path to the input GFA file representing the assembly graph.
    output_path : str
        Directory path where the output statistics file will be written.

    Returns:
    -------
    str
        Full path to the written statistics output file.

    Statistics Calculated:
    - Number of segments (nodes)
    - Number of links (edges)
    - Number of connected components
    - Number of self-loops
    - Average node degree
    - Total segment sequence length
    - Average segment length
    - N50 and L50 contiguity metrics
    - GC content across all segments
    """

    ug = UnitigGraph.from_gfa(gfa_file)

    stats = {
        "nsegments": ug.graph.vcount(),
        "nlinks": ug.graph.ecount(),
        "ncomponents": len(ug.get_connected_components()),
        "nloops": len(ug.self_loops),
        "average_node_degree": 0,
        "total_length": 0,
        "average_segment_length": 0,
        "n50": 0,
        "l50": 0,
        "gc_content": 0.0,
    }

    stats["average_node_degree"] = ug.calculate_average_node_degree()
    stats["total_length"] = ug.calculate_total_length()
    stats["average_segment_length"] = ug.calculate_average_segment_length()
    stats["n50"], stats["l50"] = ug.calculate_n50_l50()
    stats["gc_content"] = ug.get_gc_content()

    output_file = _write_stats_file(gfa_file, stats, output_path)

    # Log the statistics
    logger.info(f"Basic graph statistics for {gfa_file}:")
    logger.info(f"Number of segments: {stats['nsegments']}")
    logger.info(f"Number of links: {stats['nlinks']}")
    logger.info(f"Number of self-loops: {stats['nloops']}")
    logger.info(f"Number of connected components: {stats['ncomponents']}")
    logger.info(f"Average node degree: {stats['average_node_degree']}")
    logger.info(f"Sequence-based statistics for {gfa_file}:")
    logger.info(f"Total length of segments: {stats['total_length']} bp")
    logger.info(f"Average segment length: {stats['average_segment_length']} bp")
    logger.info(f"N50: {stats['n50']} bp")
    logger.info(f"L50: {stats['l50']} segment(s)")
    logger.info(f"GC content: {stats['gc_content']:.2%}")

    return output_file
