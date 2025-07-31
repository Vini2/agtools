#!/usr/bin/env python3

from bidict import bidict
from Bio import SeqIO
from igraph import Graph

from agtools.core.contig_graph import ContigGraph
from agtools.core.fasta_parser import FastaParser


def _get_links_megahit(gfa_file: str) -> tuple:
    """
    Parse a GFA file to extract segment sequences and connectivity (links) between segments.

    Parameters
    ----------
    gfa_file : str
        Path to the MEGAHIT-style GFA file.

    Returns
    -------
    node_count : int
        Number of unique segments.
    graph_contig_seqs : dict
        Mapping of segment ID -> sequence length in graph file.
    links : list of list
        List of 2-element lists representing linked segment IDs.
    contig_names : bidict
        Mapping of numeric node ID -> segment ID.
    """

    node_count = 0

    graph_contig_seqs = {}

    links = []

    contig_names = bidict()

    # Get links from .gfa file
    with open(gfa_file) as file:
        line = file.readline()

        while line != "":
            # Identify lines with link information
            if line.startswith("L"):
                link = []

                strings = line.split("\t")

                link1 = strings[1]
                link2 = strings[3]

                link.append(link1)
                link.append(link2)
                links.append(link)

            elif line.startswith("S"):
                strings = line.split()

                contig_names[node_count] = strings[1]

                graph_contig_seqs[strings[1]] = len(strings[2])

                node_count += 1

            line = file.readline()

    return node_count, graph_contig_seqs, links, contig_names


def _get_graph_edges_megahit(links: list, contig_names_rev: bidict) -> tuple:
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
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append((contig_names_rev[link[0]], contig_names_rev[link[1]]))
        else:
            self_loops.append(contig_names_rev[link[0]])

    return edge_list, self_loops


def get_contig_graph(gfa_file: str, contigs_file: str) -> ContigGraph:
    """
    Build a contig-level graph from a MEGAHIT GFA file and a contig FASTA file.

    Matches sequences between GFA and FASTA to map contig IDs, constructs an igraph
    representation of the graph, and packages it in a ContigGraph object.

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

    original_contig_seqs = {}
    contig_descriptions = {}

    # Get mapping of original contig identifiers with descriptions
    for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
        original_contig_seqs[record.id] = len(record.seq)
        contig_descriptions[record.id] = record.description

    # Get links and contigs of the assembly graph
    (
        node_count,
        graph_contig_seqs,
        links,
        contig_names,
    ) = _get_links_megahit(gfa_file)

    # Get list of edges and self loops
    edge_list, self_loops = _get_graph_edges_megahit(
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

    # Map original contig identifiers to contig identifiers of MEGAHIT assembly graph
    graph_to_contig_map = bidict()

    for (n, m), (n2, m2) in zip(
        graph_contig_seqs.items(), original_contig_seqs.items()
    ):
        if m == m2:
            graph_to_contig_map[n] = n2

    # Clean up temporary sequence maps
    del graph_contig_seqs
    del original_contig_seqs

    # Get parser for contigs.fasta
    parser = FastaParser(contigs_file, assembler="megahit", mapping=graph_to_contig_map)

    contig_graph = ContigGraph(
        graph=graph,
        vcount=graph.vcount(),
        ecount=graph.ecount(),
        file_path=gfa_file,
        contig_names=contig_names,
        contig_parser=parser,
        contig_descriptions=contig_descriptions,
        graph_to_contig_map=graph_to_contig_map,
        self_loops=self_loops,
    )

    return contig_graph
