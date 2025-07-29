#!/usr/bin/env python3

import tempfile
import os
import sys

from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "0.1.0"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Alpha"


def concat(graph_files: list, output_path: str) -> str:
    """
    Concatenate multiple GFA files into a single output GFA file.

    This function reads multiple GFA files, groups by tags, verifies
    the uniqueness of segments, and writes the concatenated result.

    Parameters
    ----------
    graph_files : list of str
        Paths to the GFA files to concatenate.
    output_path : str
        Directory where the concatenated GFA file will be saved.

    Returns
    -------
    str
        Path to the final concatenated GFA file.
    """

    gfa_tags = ["#", "H", "S", "L", "J", "C", "P", "W"]
    temp_files = {tag: tempfile.NamedTemporaryFile(mode="w+", delete=False) for tag in gfa_tags}
    other_lines = tempfile.NamedTemporaryFile(mode="w+", delete=False)

    output_file = f"{output_path}/concatenated_graph.gfa"

    segments = set()

    try:
        # Single pass per file: distribute lines into per-tag temp files
        for graph_file in graph_files:
            with open(graph_file, "r") as f:
                for line in f:
                    tag = line[0]

                    if tag == "S":
                        parts = line.strip().split("\t")
                        segment_id = parts[1]

                        if segment_id not in segments:
                            segments.add(segment_id)
                            
                        else:
                            logger.error("Duplicate segment IDs found in GFA files.")
                            logger.error("Please rename segment IDs and concatenate.")

                            sys.exit(1)

                    # Handle missing newline
                    line_to_write = line if line.endswith("\n") else f"{line}\n"

                    if tag in temp_files:
                        temp_files[tag].write(line_to_write)
                    else:
                        other_lines.write(line_to_write)

        # Write to concatenated output
        with open(output_file, "w") as out:
            # Write each tag group in order
            for tag in gfa_tags:
                tf = temp_files[tag]
                tf.flush()
                tf.seek(0)
                for line in tf:
                    out.write(line)

            # Write any unrecognised tags at the end
            other_lines.flush()
            other_lines.seek(0)
            for line in other_lines:
                out.write(line)

    finally:
        # Clean up temp files
        for tf in temp_files.values():
            name = tf.name
            tf.close()
            os.remove(name)
        name = other_lines.name
        other_lines.close()
        os.remove(name)
    
    return output_file
