#!/usr/bin/env python3

import io
from collections import defaultdict

from agtools.assemblers._contig_graph_base import (
    build_contig_graph as _build_contig_graph,
)
from agtools.assemblers._contig_graph_base import get_unitig_graph as _get_unitig_graph
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

    with io.open(contig_paths, mode="r", buffering=1024 * 1024) as file:
        for line in file.readlines():
            if not (line.startswith("#") or line.startswith("seq_name")):
                strings = line.strip().split()

                contig_name = strings[0]
                contig_id = len(contig_names)
                contig_name_to_id[contig_name] = contig_id
                contig_names.append(contig_name)

                path = strings[-1]
                path = path.replace("*", "")

                if path.startswith(","):
                    path = path[1:]

                if path.endswith(","):
                    path = path[:-1]

                segments = path.rstrip().split(",")

                for segment in segments:
                    if segment[0] == "-":
                        segment_contigs[segment_name_to_id[f"edge_{segment[1:]}"]].add(
                            contig_id
                        )
                    else:
                        segment_contigs[segment_name_to_id[f"edge_{segment}"]].add(
                            contig_id
                        )

    return segment_contigs, contig_names, contig_name_to_id


def get_contig_graph(
    graph_file: str, contigs_file: str, contig_paths_file: str
) -> ContigGraph:
    """
    Build a contig-level graph from an assembly GFA file and contig path mappings.

    This function parses contig metadata, links, and path structure to construct an
    undirected graph where each node represents a contig and edges represent linkages
    inferred from shared segments or GFA link data.

    Parameters
    ----------
    graph_file : str
        Path to the GFA file.
    contigs_file : str
        Path to the FASTA file with contig sequences.
    contig_paths_file : str
        Path to the file with segment paths used to build contigs.

    Returns
    -------
    ContigGraph
        An object representing the contig-level graph with node metadata.
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
