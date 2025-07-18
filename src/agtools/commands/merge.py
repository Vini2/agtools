#!/usr/bin/env python3

import sys

from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Alpha"


def _combine_gfa_files(graph_files: list) -> tuple:
    """
    Read and merge the components of multiple GFA files.

    This function parses a list of GFA files and categorizes each line by type
    (e.g., segment, link, path), ensuring segment IDs are unique across all files.
    If duplicate segment IDs are found, the program exits with an error.

    Parameters
    ----------
    graph_files : list of str
        List of paths to GFA files to be merged.

    Returns
    -------
    tuple
        Eight lists containing lines from the input GFA files, grouped as:
        (comment_lines, header_lines, segment_lines, link_lines,
         jump_lines, containment_lines, path_lines, walk_lines)
    """

    segments = {}
    comment_lines = []
    header_lines = []
    segment_lines = []
    link_lines = []
    jump_lines = []
    containment_lines = []
    path_lines = []
    walk_lines = []

    for gfa_file in graph_files:

        with open(gfa_file) as f:
            for line in f:
                if line.startswith("S"):

                    parts = line.strip().split("\t")

                    if parts[1] not in segments:
                        segments[parts[1]] = parts[2]
                        segment_lines.append(line.strip())
                    else:
                        logger.error("Duplicate segment IDs found in GFA files.")
                        logger.error("Please rename segment IDs and merge.")
                        sys.exit(1)

                elif line.startswith("L"):
                    link_lines.append(line.strip())

                elif line.startswith("J"):
                    jump_lines.append(line.strip())

                elif line.startswith("C"):
                    containment_lines.append(line.strip())

                elif line.startswith("P"):
                    path_lines.append(line.strip())

                elif line.startswith("W"):
                    walk_lines.append(line.strip())

                elif line.startswith("H"):
                    header_lines.append(line.strip())

                elif line.startswith("#"):
                    comment_lines.append(line.strip())

    return (
        comment_lines,
        header_lines,
        segment_lines,
        link_lines,
        jump_lines,
        containment_lines,
        path_lines,
        walk_lines,
    )


def _write_gfa_elements(
    comments: list,
    headers: list,
    segments: list,
    links: list,
    jumps: list,
    containments: list,
    paths: list,
    walks: list,
    output_path: str,
) -> str:
    """
    Write categorized GFA lines to a new merged GFA file.

    This function writes the collected GFA components (e.g., segments, links, paths)
    to a single output file in the correct GFA format and order.

    Parameters
    ----------
    comments : list
        Comment lines beginning with '#'.
    headers : list
        Header lines beginning with 'H'.
    segments : list
        Segment lines beginning with 'S'.
    links : list
        Link lines beginning with 'L'.
    jumps : list
        Jump lines beginning with 'J'.
    containments : list
        Containment lines beginning with 'C'.
    paths : list
        Path lines beginning with 'P'.
    walks : list
        Walk lines beginning with 'W'.
    output_path : str
        Directory where the merged GFA file will be saved.

    Returns
    -------
    str
        Path to the written merged GFA file.
    """

    output_file = f"{output_path}/merged_graph.gfa"

    with open(output_file, "w") as file_out:

        # Write comments
        for line in comments:
            file_out.write(line + "\n")

        # Write headers
        for line in headers:
            file_out.write(line + "\n")

        # Write segments
        for line in segments:
            file_out.write(line + "\n")

        # Write links
        for line in links:
            file_out.write(line + "\n")

        # Write jumps
        for line in jumps:
            file_out.write(line + "\n")

        # Write containments
        for line in containments:
            file_out.write(line + "\n")

        # Write paths
        for line in paths:
            file_out.write(line + "\n")

        # Write walks
        for line in walks:
            file_out.write(line + "\n")

    return output_file


def merge(graph_files: str, output_path: str) -> str:
    """
    Merge multiple GFA files into a single output GFA file.

    This is the main function that coordinates reading multiple GFA files,
    verifying uniqueness of segments, and writing the merged result.

    Parameters
    ----------
    graph_files : list of str
        Paths to the GFA files to merge.
    output_path : str
        Directory where the merged GFA file will be saved.

    Returns
    -------
    str
        Path to the final merged GFA file.
    """

    comments, headers, segments, links, jumps, containments, paths, walks = (
        _combine_gfa_files(graph_files)
    )
    output_file = _write_gfa_elements(
        comments,
        headers,
        segments,
        links,
        jumps,
        containments,
        paths,
        walks,
        output_path,
    )
    return output_file
