#!/usr/bin/env python3

from agtools.core.graph import UnitigGraph
from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Alpha"


def _calculate_average_node_degree(graph: UnitigGraph) -> int:
    """
    Calculate the average node degree of the graph.

    Parameters
    ----------
    graph : UnitigGraph
        The unitig graph object containing the assembly graph.

    Returns
    -------
    int
        Average node degree of the graph.
    """

    if graph.graph.vcount() == 0:
        return 0
    return int(sum(graph.graph.degree()) / graph.graph.vcount())


def _calculate_total_length(segment_lengths: dict) -> int:
    """
    Calculate the total length of all segments in the graph.

    Parameters
    ----------
    segment_lengths : dict
        Dictionary mapping segment IDs to their lengths.

    Returns
    -------
    int
        Total length of all segments.
    """
    return sum(segment_lengths.values())


def _calculate_average_segment_length(segment_lengths: dict) -> int:
    """
    Calculate the average segment length.

    Parameters
    ----------
    segment_lengths : dict
        Dictionary mapping segment IDs to their lengths.

    Returns
    -------
    int
        Average segment length.
    """
    if not segment_lengths:
        return 0
    return int(sum(segment_lengths.values()) / len(segment_lengths))


def _calculate_n50_l50(lengths: list[int]) -> tuple[int, int]:
    """
    Calculate N50 and L50 from a list of segment lengths.

    Parameters
    ----------
    lengths : list of int
        List of segment lengths.

    Returns
    -------
    tuple of (int, int)
        A tuple containing:
        - N50 : int
            The length N such that 50% of the total length is contained in segments of length ≥ N.
        - L50 : int
            The minimum number of segments whose summed length ≥ 50% of the total.
    """
    if not lengths:
        return (0, 0)

    sorted_lengths = sorted(lengths, reverse=True)
    total_length = sum(sorted_lengths)
    cum_sum = 0

    for i, length in enumerate(sorted_lengths):
        cum_sum += length
        if cum_sum >= total_length / 2:
            return length, i + 1


def _get_gc_content(sequences: list, total_length: int) -> float:
    """
    Calculate the GC content of sequences.

    Parameters
    ----------
    sequence : list of str
        A list of nucleotide sequences (A, T, G, C).

    Returns
    -------
    float
        GC content as a percentage of total base pairs.
    """
    if not sequences:
        return 0.0
    elif total_length == 0:
        return 0.0

    gc_count = sum(seq.count("G") + seq.count("C") for seq in sequences)
    return gc_count / total_length


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
        "ncomponents": len(ug.graph.components()),
        "nloops": len(ug.self_loops),
        "average_node_degree": 0,
        "total_length": 0,
        "average_segment_length": 0,
        "n50": 0,
        "l50": 0,
        "gc_content": 0.0,
    }

    stats["average_node_degree"] = _calculate_average_node_degree(ug)
    stats["total_length"] = _calculate_total_length(ug.segment_lengths)
    stats["average_segment_length"] = _calculate_average_segment_length(
        ug.segment_lengths
    )
    stats["n50"], stats["l50"] = _calculate_n50_l50(ug.segment_lengths.values())
    stats["gc_content"] = _get_gc_content(
        ug.segment_sequences.values(), stats["total_length"]
    )

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
