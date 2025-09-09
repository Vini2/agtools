#!/usr/bin/env python3

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "1.0.1"
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
    links = []

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
                parts = line.strip().split("\t")[1].split(" ")
                seq1_name = parts[0]
                seq2_name = parts[1]
                seq1_overlap = int(parts[3]) - int(parts[2])
                seq2_overlap = int(parts[6]) - int(parts[5])
                seq2_orient = int(parts[8])

                if seq1_overlap == seq2_overlap:

                    # seq2 is reversed with respect to seq1
                    if seq2_orient == 1:
                        links.append([seq1_name, "+", seq2_name, "-", seq1_overlap])
                    elif seq2_orient == 0:
                        links.append([seq1_name, "+", seq2_name, "+", seq1_overlap])

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
        Directory path where the output GFA file will be saved.

    Returns
    -------
    str
        Path to the generated GFA file.
    """

    output_file = f"{output_path}/converted_graph.gfa"

    with open(output_file, "w") as gfa_file:

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
        Directory path where the output GFA file will be saved.

    Returns
    -------
    str
        Path to the converted GFA file.
    """

    segments, links = _get_segments_and_links(asqg_file)

    output_file = _write_gfa(segments, links, output_path)

    return output_file
