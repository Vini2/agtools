#!/usr/bin/env python3

from agtools.commands._output import prepare_output_file
from agtools.core.gfa_filter import write_filtered_gfa
from agtools.core.unitig_graph import UnitigGraph
from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "1.0.2"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def _write_filtered_graph(
    segments_to_remove: set, gfa_file: str, output_path: str
) -> str:
    """
    Write a filtered GFA file by excluding lines that involve specified segments.

    This function processes a GFA file line by line, excluding:
    - `S` lines (segments) whose ID is in `segments_to_remove`
    - `L`, `J`, and `C` lines (links, joins, containments) involving any removed segments
    - `P` and `W` lines (paths, walks) that contain any removed segment IDs

    All other lines (including headers or comments) are preserved.

    Parameters
    ----------
    segments_to_remove : set[str]
        Segment IDs that should be removed from the graph.
    gfa_file : str
        Path to the input GFA file.
    output_path : str
        Path where the filtered GFA file will be written.

    Returns
    -------
    str
        Full path to the written filtered GFA file.
    """

    output_file = prepare_output_file(output_path)

    def keep_segment(seg_id: str) -> bool:
        return seg_id not in segments_to_remove

    write_filtered_gfa(gfa_file, output_file, keep_segment)

    return output_file


def filter(gfa_file: str, min_length: int, output_path: str) -> str:
    """
    Filter segments from a GFA file based on minimum sequence length.

    This function loads a GFA-formatted assembly graph using UnitigGraph,
    identifies all segments shorter than the specified `min_length`, and
    writes a new filtered GFA file excluding those segments and all lines
    (links, paths, etc.) involving them.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.
    min_length : int
        Minimum length threshold for retaining segments.
    output_path : str
        Path where the filtered GFA file will be saved.

    Returns
    -------
    str
        Full path to the filtered GFA file.
    """

    ug = UnitigGraph.from_gfa(gfa_file)

    segments_to_remove = {
        seg_id
        for seg_id, seq_length in ug.segment_lengths.items()
        if seq_length < min_length
    }
    logger.info(
        f"Identified {len(segments_to_remove)} segments shorter than {min_length} bp to remove"
    )

    output_file = _write_filtered_graph(segments_to_remove, gfa_file, output_path)

    return output_file


# TODO: filter isolated segments
