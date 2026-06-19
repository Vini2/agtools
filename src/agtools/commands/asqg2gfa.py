#!/usr/bin/env python3

from agtools import __version__
from agtools.commands._format_checks import validate_asqg_input
from agtools.commands._output import open_output_file
from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = __version__
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def _get_segments_and_links(asqg_file: str) -> tuple:
    """
    Parses an ASQG (Assembly String Graph) file to extract segments and links.

    Parameters
    ----------
    asqg_file : str
        Path to the input ASQG file.

    Returns
    -------
    tuple
        A tuple containing:
        - segments (dict): Mapping of segment ID to sequence string.
        - links (list): List of links, each represented as
          [from_segment, from_orientation, to_segment, to_orientation, overlap_length].

    References
    ----------
    ASQG Format
    [https://github.com/jts/sga/wiki/ASQG-Format](https://github.com/jts/sga/wiki/ASQG-Format)
    """

    segments = {}
    raw_edges = []

    # Get contig connections from .asqg file
    with open(asqg_file) as file:
        for line in file.readlines():

            # Count the number of contigs
            if line.startswith("VT"):
                parts = line.strip().split("\t")
                contig_name = parts[1]
                contig_seq = parts[2]
                segments[contig_name] = contig_seq

            # Identify lines with link information
            elif line.startswith("ED"):
                fields = line.strip().split("\t")
                if len(fields) < 2:
                    message = f"Malformed ED line: {line.strip()}"
                    logger.error(message)
                    raise ValueError(message)

                parts = fields[1].split(" ")
                if len(parts) < 10:
                    message = f"Malformed ED line: {line.strip()}"
                    logger.error(message)
                    raise ValueError(message)

                seq1_name = parts[0]
                seq2_name = parts[1]

                try:
                    seq1_start = int(parts[2])
                    seq1_end = int(parts[3])
                    seq1_length = int(parts[4])
                    seq2_start = int(parts[5])
                    seq2_end = int(parts[6])
                    seq2_length = int(parts[7])
                    seq2_orient = int(parts[8])
                    overlap_dif = int(parts[9])
                except ValueError as e:
                    message = f"Malformed ED line: {line.strip()}"
                    logger.error(message)
                    raise ValueError(message) from e

                raw_edges.append(
                    [
                        seq1_name,
                        seq2_name,
                        seq1_start,
                        seq1_end,
                        seq1_length,
                        seq2_start,
                        seq2_end,
                        seq2_length,
                        seq2_orient,
                        line.strip(),
                    ]
                )

    links = []
    for (
        seq1_name,
        seq2_name,
        seq1_start,
        seq1_end,
        seq1_length,
        seq2_start,
        seq2_end,
        seq2_length,
        seq2_orient,
        raw_line,
    ) in raw_edges:
        if seq1_name not in segments or seq2_name not in segments:
            message = f"Malformed ED line: {raw_line}"
            logger.error(message)
            raise ValueError(message)

        seq1_overlap = seq1_end - seq1_start + 1
        seq2_overlap = seq2_end - seq2_start + 1

        if (
            seq1_start < 0
            or seq1_end < seq1_start
            or seq1_end >= seq1_length
            or seq2_start < 0
            or seq2_end < seq2_start
            or seq2_end >= seq2_length
            or seq1_overlap != seq2_overlap
            or seq1_length != len(segments[seq1_name])
            or seq2_length != len(segments[seq2_name])
        ):
            message = f"Malformed ED line: {raw_line}"
            logger.error(message)
            raise ValueError(message)

        overlap = seq1_overlap
        seq1_is_suffix = (
            seq1_start == seq1_length - overlap and seq1_end == seq1_length - 1
        )
        seq1_is_prefix = seq1_start == 0 and seq1_end == overlap - 1

        if seq1_is_suffix:
            seq1_orient = "+"
        elif seq1_is_prefix:
            seq1_orient = "-"
        else:
            message = f"Malformed ED line: {raw_line}"
            logger.error(message)
            raise ValueError(message)

        if seq2_orient == 0:
            if not (seq2_start == 0 and seq2_end == overlap - 1):
                message = f"Malformed ED line: {raw_line}"
                logger.error(message)
                raise ValueError(message)
            seq2_gfa_orient = "+"
        elif seq2_orient == 1:
            if not (
                seq2_start == seq2_length - overlap and seq2_end == seq2_length - 1
            ):
                message = f"Malformed ED line: {raw_line}"
                logger.error(message)
                raise ValueError(message)
            seq2_gfa_orient = "-"
        else:
            message = f"Malformed ED line: {raw_line}"
            logger.error(message)
            raise ValueError(message)

        links.append([seq1_name, seq1_orient, seq2_name, seq2_gfa_orient, overlap])

    return segments, links


def _write_gfa(segments, links, output_path):
    """
    Writes segments and links to a GFA (Graphical Fragment Assembly) file.

    Parameters
    ----------
    segments : dict
        Dictionary of segment IDs to sequences.
    links : list
        List of link definitions in the form
        [from_segment, from_orientation, to_segment, to_orientation, overlap_length].
    output_path : str
        Path where the output GFA file will be saved.

    Returns
    -------
    str
        Path to the generated GFA file.
    """

    with open_output_file(output_path) as (output_file, gfa_file):

        # Write segments
        for seg_id, seq in segments.items():
            gfa_file.write(f"S\t{seg_id}\t{seq}\n")

        # Write links
        for link in links:
            from_seg, from_orient, to_seg, to_orient, overlap = link
            orient_str = "+" if from_orient == "+" else "-"
            gfa_file.write(
                f"L\t{from_seg}\t{orient_str}\t{to_seg}\t{to_orient}\t{overlap}M\n"
            )

    return output_file


def asqg2gfa(asqg_file, output_path):
    """
    Converts an ASQG file to a GFA file.

    This function parses segment and link data from an ASQG file and writes them
    into a GFA-format file for downstream graph analysis or visualization.

    Parameters
    ----------
    asqg_file : str
        Path to the input ASQG file.
    output_path : str
        Path where the output GFA file will be saved.

    Returns
    -------
    str
        Path to the converted GFA file.
    """

    validate_asqg_input(asqg_file, "asqg2gfa")
    segments, links = _get_segments_and_links(asqg_file)

    output_file = _write_gfa(segments, links, output_path)

    return output_file
