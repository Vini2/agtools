#!/usr/bin/env python3

from agtools import __version__
from agtools.commands._output import prepare_output_file
from agtools.core.fasta_parser import FastaParser
from agtools.core.gfa_filter import write_filtered_gfa
from agtools.core.unitig_graph import UnitigGraph

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = __version__
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def _write_filtered_graph(
    segments_to_remove: set, parser: FastaParser, gfa_file: str, output_path: str
) -> str:
    """
    Write a cleaned GFA file by excluding lines that involve specified segments.

    This function processes a GFA file line by line, excluding:
    - `S` lines (segments) whose ID is in `segments_to_remove`
    - `L`, `J`, and `C` lines (links, joins, containments) involving any removed segments
    - `P` and `W` lines (paths, walks) that contain any removed segment IDs

    All other lines (including headers or comments) are preserved.

    Parameters
    ----------
    segments_to_remove : set[str]
        Segment IDs that should be removed from the graph.
    parser : FastaParser
        FASTA parser.
    gfa_file : str
        Path to the input GFA file.
    output_path : str
        Path where the cleaned GFA file will be written.

    Returns
    -------
    str
        Full path to the written cleaned GFA file.
    """

    output_file = prepare_output_file(output_path)

    def keep_segment(seg_id: str) -> bool:
        return seg_id not in segments_to_remove

    def transform_segment(parts: list[str]) -> list[str]:
        if len(parts) > 2 and parts[2] in ("", "*"):
            updated_parts = parts.copy()
            updated_parts[2] = str(parser.get_sequence(parts[1]))
            return updated_parts
        return parts

    write_filtered_gfa(gfa_file, output_file, keep_segment, transform_segment)

    return output_file


def clean(gfa_file: str, fasta: str, assembler: str, output_path: str) -> str:
    """
    Clean a GFA file based on segments in a FASTA file.

    This function adds the FASTA sequence to the GFA file if
    missing, removes segments if not present in the FASTA file
    and removes any links, paths, walks, junctions or
    containments containing missing segments.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.
    fasta : int
        Path to the FASTA file.
    assembler : str
        Assembler used to get the assembly
    output_path : str
        Path where the filtered GFA file will be saved.

    Returns
    -------
    str
        Full path to the cleaned GFA file.
    """

    ug = UnitigGraph.from_gfa(gfa_file)

    # Get parser for fasta file
    parser = FastaParser(fasta, assembler=assembler)

    segments_to_remove = {
        segment for segment in ug.segment_names if segment not in parser.index
    }

    output_file = _write_filtered_graph(
        segments_to_remove, parser, gfa_file, output_path
    )

    return output_file
