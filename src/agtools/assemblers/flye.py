#!/usr/bin/env python3

import io
from collections import defaultdict

from igraph import Graph

from agtools.core.contig_graph import ContigGraph
from agtools.core.fasta_parser import FastaParser
from agtools.core.unitig_graph import UnitigGraph


def _get_segments(graph_file: str) -> dict:
    """
    Parse a GFA file to extract segment names and their corresponding IDs.

    Parameters
    ----------
    graph_file : str
        Path to the GFA file.

    Returns
    -------
    segment_name_to_id : dict[str, int]
        Mapping from segment name to its internal ID.
    """

    segment_name_to_id = dict()
    segment_names = list()

    with io.open(graph_file, mode="r", buffering=1024 * 1024) as f:
        while True:
            line = f.readline()
            if not line:
                break

            tag = line[0]

            if not line:
                continue

            if tag == "S":  # Segment line
                parts = line.rstrip().split("\t")
                seg_name = parts[1]
                seg_id = len(segment_names)
                segment_name_to_id[seg_name] = seg_id
                segment_names.append(seg_name)

    return segment_name_to_id


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
    contig_id = -1

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


def _get_graph_edges(
    graph_file: str, segment_contigs: dict, segment_name_to_id: dict
) -> tuple:
    """
    Construct edges between contigs based on shared segment links in the GFA file.

    Parameters
    ----------
    graph_file : str
        Path to the GFA file.
    segment_contigs : dict[str, set[str]]
        Mapping from segment ID to contigs containing them.
    segment_name_to_id : dict[str, int]
        Mapping from segment name to its internal ID.

    Returns
    -------
    tuple
    list : tuple[int, int]
        List of edges as (source_node_id, target_node_id).
    self_loops : list
        List of self loops
    lcount : int
        The number of links (lines starting with tag "L") in the graph
    """

    lcount = 0
    self_loops = set()
    edge_list = set()

    # Get links from assembly_graph_with_scaffolds.gfa
    with io.open(graph_file, mode="r", buffering=1024 * 1024) as file:
        line = file.readline()

        while line != "":
            # Identify lines with link information
            if "L" in line:
                lcount += 1
                strings = line.split("\t")
                source = segment_name_to_id[strings[1]]
                target = segment_name_to_id[strings[3]]

                source_contigs = None
                target_contigs = None

                if source in segment_contigs:
                    source_contigs = segment_contigs[source]

                if target in segment_contigs:
                    target_contigs = segment_contigs[target]

                if source_contigs and target_contigs:
                    for source_contig in source_contigs:
                        for target_contig in target_contigs:
                            if (
                                source_contig != target_contig
                                and (source_contig, target_contig) not in edge_list
                            ):
                                edge_list.add((source_contig, target_contig))
                            else:
                                self_loops.add(source_contig)

            line = file.readline()

    return list(edge_list), list(self_loops), lcount


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

    # Get segment names and their IDs from the GFA file
    segment_name_to_id = _get_segments(graph_file)

    # Get paths, segments, links and contigs of the assembly graph
    (segment_contigs, contig_names, contig_name_to_id) = (
        _get_segment_paths_and_contig_mapping(contig_paths_file, segment_name_to_id)
    )

    node_count = len(contig_names)

    # Create graph
    graph = Graph()

    # Add vertices
    graph.add_vertices(node_count)

    # Name vertices with contig identifiers
    graph.vs["label"] = contig_names

    # Get list of edges
    edge_list, self_loops, lcount = _get_graph_edges(
        graph_file=graph_file,
        segment_contigs=segment_contigs,
        segment_name_to_id=segment_name_to_id,
    )

    # Add edges to the graph
    graph.add_edges(edge_list)

    # Simplify the graph
    graph.simplify(multiple=True, loops=False, combine_edges=None)

    # Get parser for contigs.fasta
    parser = FastaParser(contigs_file)

    # Create ContigGraph object
    contig_graph = ContigGraph(
        graph=graph,
        vcount=graph.vcount(),
        lcount=lcount,
        ecount=graph.ecount(),
        file_path=graph_file,
        contig_names=contig_names,
        contig_name_to_id=contig_name_to_id,
        contig_parser=parser,
        contig_descriptions=None,
        graph_to_contig_map=None,
        self_loops=self_loops,
    )

    return contig_graph


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

    ug = UnitigGraph.from_gfa(graph_file)
    return ug
