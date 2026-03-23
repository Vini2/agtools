#!/usr/bin/env python3

import io
import re
from collections import defaultdict

from agtools.assemblers._contig_graph_base import (
    build_contig_graph as _build_contig_graph,
    get_unitig_graph as _get_unitig_graph,
)
from agtools.core.contig_graph import ContigGraph
from agtools.core.unitig_graph import UnitigGraph

def _get_segment_paths_and_contig_mapping(
    contig_paths: str, segment_name_to_id: dict
) -> tuple:
    """
    Parse a contig paths file and extract segment-contig relationships.

    Parameters
    ----------
    contig_paths : str
        Path to the contig paths file (e.g. contigs.paths of scaffolds.paths).
    segment_name_to_id : dict[str, int]
        Mapping from segment name to its internal ID.

    Returns
    -------
    tuple
    segment_contigs : dict[str, set[str]]
        Mapping from segment ID to the set of contig numbers it appears in.
    contig_names : list
        List of contig names.
    contig_name_to_id : dict[str, int]
        Mapping from contig name to its internal ID.
    """

    contig_names = []
    contig_name_to_id = dict()

    segment_contigs = defaultdict(set)

    current_contig_num = -1

    with io.open(contig_paths, mode="r", buffering=1024 * 1024) as file:
        name = file.readline().rstrip()
        path = file.readline().rstrip("\n")

        while name != "" and path != "":
            while ";" in path:
                path = path[:-1] + "," + file.readline().rstrip("\n")

            start = "NODE_"
            end = "_length_"
            contig_num = int(re.search("%s(.*)%s" % (start, end), name).group(1))

            segments = path.rstrip().split(",")

            if current_contig_num != contig_num:
                segment_ids = [segment_name_to_id[seg[:-1]] for seg in segments]

                contig_id = len(contig_names)
                contig_name_to_id[name] = contig_id
                contig_names.append(name)
                current_contig_num = contig_num

                for segment_id in segment_ids:
                    segment_contigs[segment_id].add(contig_id)

            name = file.readline().rstrip()
            path = file.readline().rstrip()

    return segment_contigs, contig_names, contig_name_to_id


def get_contig_graph(
    graph_file: str, contigs_file: str, contig_paths_file: str
) -> ContigGraph:
    """
    Build a contig-level graph from a GFA file and a contig paths mapping file.

    Parameters
    ----------
    graph_file : str
        Path to the GFA file.
    contigs_file : str
        Path to the FASTA file with contig sequences.
    contig_paths_file : str
        Path to the contigs.paths or scaffolds.paths file.

    Returns
    -------
    ContigGraph
        Parsed contig graph object.
    """

    return _build_contig_graph(
        graph_file=graph_file,
        contigs_file=contigs_file,
        contig_paths_file=contig_paths_file,
        segment_path_parser=_get_segment_paths_and_contig_mapping,
    )


def get_unitig_graph(graph_file: str) -> UnitigGraph:
    """
    Build a unitig-level assembly graph from a GFA file.

    Parameters
    ----------
    graph_file : str
        Path to the GFA file.

    Returns
    -------
    UnitigGraph
        Parsed unitig graph object.
    """

    return _get_unitig_graph(graph_file)
