#!/usr/bin/env python3

import sys

from agtools.log_config import logger


def _combine_gfa_files(graph_files):

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

    return comment_lines, header_lines, segment_lines, link_lines, jump_lines, containment_lines, path_lines, walk_lines


def _write_gfa_elements(comments, headers, segments, links, jumps, containments, paths, walks, output_path):
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


def merge(graph_files, output_path):
    comments, headers, segments, links, jumps, containments, paths, walks = _combine_gfa_files(graph_files)
    output_file = _write_gfa_elements(comments, headers, segments, links, jumps, containments, paths, walks, output_path)
    return output_file
