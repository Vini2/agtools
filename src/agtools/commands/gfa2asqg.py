#!/usr/bin/env python3

import re

from agtools import __version__
from agtools.commands._format_checks import validate_gfa_input
from agtools.commands._output import open_output_file
from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"

__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


_OVERLAP_PATTERN = re.compile(r"^(\d+)M$")


def _raise_logged_value_error(message: str) -> None:
    logger.error(message)
    raise ValueError(message)


def _get_segment_length(parts: list[str], line: str) -> int:
    """Get the segment length from the sequence field or an ``LN`` tag."""

    sequence = parts[2]
    if sequence != "*":
        return len(sequence)

    for field in parts[3:]:
        if field.startswith("LN:i:"):
            try:
                return int(field.split(":")[-1])
            except ValueError as exc:
                message = f"Malformed S line: {line.strip()}"
                logger.error(message)
                raise ValueError(message) from exc

    _raise_logged_value_error(
        f"Cannot convert segment {parts[1]} without sequence data or an LN tag."
    )


def _get_segments_and_edges(gfa_file: str) -> tuple[dict[str, str], list[list]]:
    """
    Parse a GFA file and extract segments and ASQG-compatible edges.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.

    Returns
    -------
    tuple
        A tuple containing:
        - segments (dict): Mapping of segment ID to sequence string.
        - edges (list): List of ASQG edge fields in the order expected by ``ED`` records.
    """

    segments = {}
    segment_lengths = {}
    raw_links = []

    with open(gfa_file) as file:
        for line in file:
            if line.startswith("S\t"):
                parts = line.strip().split("\t")
                if len(parts) < 3:
                    _raise_logged_value_error(f"Malformed S line: {line.strip()}")

                segment_id = parts[1]
                segments[segment_id] = parts[2]
                segment_lengths[segment_id] = _get_segment_length(parts, line)

            elif line.startswith("L\t"):
                parts = line.strip().split("\t")
                if len(parts) < 6:
                    _raise_logged_value_error(f"Malformed L line: {line.strip()}")

                from_segment, from_orient = parts[1], parts[2]
                to_segment, to_orient = parts[3], parts[4]
                overlap_field = parts[5]

                if from_orient not in {"+", "-"} or to_orient not in {"+", "-"}:
                    _raise_logged_value_error(f"Malformed L line: {line.strip()}")

                overlap_match = _OVERLAP_PATTERN.fullmatch(overlap_field)
                if overlap_match is None:
                    _raise_logged_value_error(
                        f"Unsupported L overlap field in line: {line.strip()}"
                    )

                raw_links.append(
                    (
                        from_segment,
                        from_orient,
                        to_segment,
                        to_orient,
                        int(overlap_match.group(1)),
                    )
                )

    edges = []
    for from_segment, from_orient, to_segment, to_orient, overlap in raw_links:
        if from_segment not in segment_lengths or to_segment not in segment_lengths:
            _raise_logged_value_error(
                f"Link references undefined segment(s): {from_segment}, {to_segment}"
            )

        from_length = segment_lengths[from_segment]
        to_length = segment_lengths[to_segment]

        if overlap > from_length or overlap > to_length:
            _raise_logged_value_error(
                f"Link overlap {overlap} exceeds segment length in "
                f"{from_segment}->{to_segment}."
            )

        from_start = from_length - overlap if from_orient == "+" else 0
        from_end = from_length - 1 if from_orient == "+" else overlap - 1
        to_start = 0 if to_orient == "+" else to_length - overlap
        to_end = overlap - 1 if to_orient == "+" else to_length - 1
        to_asqg_orient = 1 if to_orient == "-" else 0

        edges.append(
            [
                from_segment,
                to_segment,
                from_start,
                from_end,
                from_length,
                to_start,
                to_end,
                to_length,
                to_asqg_orient,
                0,
            ]
        )

    return segments, edges


def _write_asqg(segments: dict[str, str], edges: list[list], output_path: str) -> str:
    """
    Write segments and edges to an ASQG file.

    Parameters
    ----------
    segments : dict[str, str]
        Mapping of segment IDs to sequence strings.
    edges : list[list]
        Edge records formatted for ASQG ``ED`` lines.
    output_path : str
        Path where the output ASQG file should be saved.

    Returns
    -------
    str
        Path to the generated ASQG file.
    """

    min_overlap = min((edge[3] - edge[2] + 1 for edge in edges), default=0)

    with open_output_file(output_path) as (output_file, asqg_file):
        asqg_file.write(f"HT\tVN:i:1\tER:f:0\tOL:i:{min_overlap}\tCN:i:0\tTE:i:0\n")

        for segment_id, sequence in segments.items():
            asqg_file.write(f"VT\t{segment_id}\t{sequence}\n")

        for edge in edges:
            asqg_file.write(
                "ED\t"
                f"{edge[0]} {edge[1]} {edge[2]} {edge[3]} {edge[4]} "
                f"{edge[5]} {edge[6]} {edge[7]} {edge[8]} {edge[9]}\n"
            )

    return output_file


def gfa2asqg(gfa_file: str, output_path: str) -> str:
    """
    Convert a GFA file to an ASQG file.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.
    output_path : str
        Path where the output ASQG file will be saved.

    Returns
    -------
    str
        Path to the converted ASQG file.
    """

    validate_gfa_input(gfa_file, "gfa2asqg")
    segments, edges = _get_segments_and_edges(gfa_file)

    output_file = _write_asqg(segments, edges, output_path)

    return output_file
