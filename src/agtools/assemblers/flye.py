#!/usr/bin/env python3

from collections import defaultdict

from bidict import bidict
from igraph import Graph

from agtools.core.contig_graph import ContigGraph
from agtools.core.fasta_parser import FastaParser
from agtools.core.unitig_graph import UnitigGraph


def _get_links(contig_paths_file: str) -> tuple:
    """
    Parse contig paths file to extract paths and segment-contig mappings.

    Parameters
    ----------
    contig_paths_file : str
        Path to the file containing contig path information.

    Returns
    -------
    tuple
        A tuple of:
        - contig_names : dict
            Bidict mapping contig indices to contig names.
        - paths : dict
            Mapping from contig index to list of segment identifiers.
        - segment_contigs : dict
            Mapping from segment ID to a set of contig indices it belongs to.
    """

    contig_names = bidict()
    contig_num = 0

    paths = {}
    segment_contigs = {}

    with open(contig_paths_file, "r") as file:
        for line in file.readlines():
            if not (line.startswith("#") or line.startswith("seq_name")):
                strings = line.strip().split()

                contig_name = strings[0]

                path = strings[-1]
                path = path.replace("*", "")

                if path.startswith(","):
                    path = path[1:]

                if path.endswith(","):
                    path = path[:-1]

                segments = path.rstrip().split(",")

                contig_names[contig_num] = contig_name

                if contig_num not in paths:
                    paths[contig_num] = segments

                for segment in segments:
                    if segment not in segment_contigs:
                        segment_contigs[segment] = set([contig_num])
                    else:
                        segment_contigs[segment].add(contig_num)

                contig_num += 1

    return contig_names, paths, segment_contigs


def _get_graph_edges(graph_file: str, paths: dict, segment_contigs: dict) -> list:
    """
    Generate edges for a contig-level graph based on GFA link information and segment-to-contig mappings.

    Parameters
    ----------
    graph_file : str
        Path to the GFA file containing link (`L`) lines.
    paths : dict
        Mapping from contig index to list of segment IDs in the path.
    segment_contigs : dict
        Mapping from segment ID to set of contig indices that contain it.

    Returns
    -------
    list of tuple
        List of edges as (source_index, target_index) pairs.
    """

    links_map = defaultdict(set)

    # Get links from assembly_graph_with_scaffolds.gfa
    with open(graph_file) as file:
        line = file.readline()

        while line != "":
            # Identify lines with link information
            if "L" in line:
                strings = line.split("\t")

                f1, f2 = "", ""

                if strings[2] == "+":
                    f1 = strings[1][5:]
                if strings[2] == "-":
                    f1 = "-" + strings[1][5:]
                if strings[4] == "+":
                    f2 = strings[3][5:]
                if strings[4] == "-":
                    f2 = "-" + strings[3][5:]

                links_map[f1].add(f2)
                links_map[f2].add(f1)

            line = file.readline()

    # Create list of edges
    edge_list = []

    for i in paths:
        segments = paths[i]

        new_links = []

        for segment in segments:
            my_segment = segment
            my_segment_num = ""

            my_segment_rev = ""

            if my_segment.startswith("-"):
                my_segment_rev = my_segment[1:]
                my_segment_num = my_segment[1:]
            else:
                my_segment_rev = "-" + my_segment
                my_segment_num = my_segment

            if my_segment in links_map:
                new_links.extend(list(links_map[my_segment]))

            if my_segment_rev in links_map:
                new_links.extend(list(links_map[my_segment_rev]))

            if my_segment in segment_contigs:
                for contig in segment_contigs[my_segment]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

            if my_segment_rev in segment_contigs:
                for contig in segment_contigs[my_segment_rev]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

            if my_segment_num in segment_contigs:
                for contig in segment_contigs[my_segment_num]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

        for new_link in new_links:
            if new_link in segment_contigs:
                for contig in segment_contigs[new_link]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

            if new_link.startswith("-"):
                if new_link[1:] in segment_contigs:
                    for contig in segment_contigs[new_link[1:]]:
                        if i != contig:
                            # Add edge to list of edges
                            edge_list.append((i, contig))

    return edge_list


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

    # Get contigs map, links and contigs of the assembly graph
    (
        contig_names,
        paths,
        segment_contigs,
    ) = _get_links(contig_paths_file)
    node_count = len(contig_names)

    # Create graph
    graph = Graph()

    # Add vertices
    graph.add_vertices(node_count)

    # Name vertices with contig identifiers
    for i in range(node_count):
        graph.vs[i]["id"] = i
        graph.vs[i]["label"] = contig_names[i]

    edge_list = _get_graph_edges(
        graph_file=graph_file,
        paths=paths,
        segment_contigs=segment_contigs,
    )

    # Add edges to the graph
    graph.add_edges(edge_list)

    # Simplify the graph
    graph.simplify(multiple=True, loops=False, combine_edges=None)

    # Get parser for contigs.fasta
    parser = FastaParser(contigs_file)

    contig_graph = ContigGraph(
        graph=graph,
        vcount=graph.vcount(),
        ecount=graph.ecount(),
        file_path=graph_file,
        contig_names=contig_names,
        contig_parser=parser,
    )

    return contig_graph


def get_unitig_graph(graph_file) -> UnitigGraph:
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
