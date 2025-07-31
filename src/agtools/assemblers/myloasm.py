#!/usr/bin/env python3

from bidict import bidict
from igraph import Graph

from agtools.core.contig_graph import ContigGraph
from agtools.core.fasta_parser import FastaParser


def _get_links_myloasm(gfa_file: str, contig_index: dict) -> tuple:
    """
    Parse a GFA file to extract contig information and connectivity
    information (links) between contigs.

    Parameters
    ----------
    gfa_file : str
        Path to the myloasm-style GFA file.

    Returns
    -------
    node_count : int
        Number of unique segments.
    links : list of list
        List of 2-element lists representing linked segment IDs.
    contig_names : bidict
        Mapping of numeric node ID -> contig ID.
    """
    node_count = 0

    links = []

    contig_names = bidict()

    # Get links from .gfa file
    with open(gfa_file) as file:
        for line in file.readlines():
            # Identify lines with link information
            if line.startswith("L"):
                link = []

                strings = line.strip().split("\t")

                link1 = strings[1]
                link1_orient = strings[2]
                link2 = strings[3]
                link2_orient = strings[4]

                if link1 in contig_index and link2 in contig_index:

                    link.append(link1)
                    link.append(link1_orient)
                    link.append(link2)
                    link.append(link2_orient)

                    links.append(link)

            # Identify lines with contig information
            elif line.startswith("S"):
                strings = line.strip().split("\t")

                if strings[1] in contig_index:
                    contig_names[node_count] = strings[1]
                    node_count += 1

    return node_count, links, contig_names


def _get_graph_edges_myloasm(links: list, contig_names_rev: bidict) -> tuple:
    """
    Convert a list of segment links into igraph-compatible edges.

    Parameters
    ----------
    links : list of list
        Pairs of linked segment IDs.
    contig_names_rev : bidict
        Mapping of segment ID -> numeric node ID.

    Returns
    -------
    edge_list : list of tuple
        List of edges as tuples of node IDs.
    self_loops : list of int
        List of node IDs that form self-loops.
    """

    edge_list = []
    self_loops = []

    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[2]:
            # Add edge to list of edges
            edge_list.append((contig_names_rev[link[0]], contig_names_rev[link[2]]))
        else:
            self_loops.append(contig_names_rev[link[0]])

    return edge_list, self_loops


def get_contig_graph(gfa_file: str, contigs_file: str) -> ContigGraph:
    """
    Build a contig-level graph from a myloasm GFA file and a contig FASTA file.

    Parameters
    ----------
    gfa_file : str
        Path to the GFA file.
    contigs_file : str
        Path to the contigs FASTA file.

    Returns
    -------
    ContigGraph
        Parsed contig graph object.
    """

    # Get parser for contigs.fasta
    parser = FastaParser(contigs_file, assembler="myloasm")

    # Get links and contigs of the assembly graph
    node_count, links, contig_names = _get_links_myloasm(gfa_file, parser.index)

    # Get list of edges and self loops
    edge_list, self_loops = _get_graph_edges_myloasm(
        links=links, contig_names_rev=contig_names.inverse
    )

    # Create graph
    graph = Graph()

    # Add vertices
    graph.add_vertices(node_count)

    # Name vertices with contig identifiers
    for i in range(node_count):
        graph.vs[i]["id"] = i
        graph.vs[i]["label"] = contig_names[i]

    # Add edges to the graph
    graph.add_edges(edge_list)

    # Simplify the graph
    graph.simplify(multiple=True, loops=False, combine_edges=None)

    contig_graph = ContigGraph(
        graph=graph,
        vcount=graph.vcount(),
        ecount=graph.ecount(),
        file_path=gfa_file,
        contig_names=contig_names,
        contig_parser=parser,
        contig_descriptions=None,
        graph_to_contig_map=None,
        self_loops=self_loops,
    )

    return contig_graph
