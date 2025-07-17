#!/usr/bin/env python3

from bidict import bidict
from Bio import SeqIO
from igraph import Graph

from agtools.core.graph import ContigGraph


def _get_links_megahit(gfa_file):
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

                graph_contig_seqs[strings[1]] = strings[2]

                node_count += 1

            line = file.readline()

    return node_count, graph_contig_seqs, links, contig_names


def _get_graph_edges_megahit(links, contig_names_rev):
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


def get_contig_graph(gfa_file, contigs_file):

    graph_contig_seqs = {}
    contig_descriptions = {}
    contig_sequences = {}

    # Get mapping of original contig identifiers with descriptions
    for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
        contig_sequences[record.id] = record.seq
        graph_contig_seqs[record.id] = str(record.seq)
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

    for (n, m), (n2, m2) in zip(graph_contig_seqs.items(), graph_contig_seqs.items()):
        if m == m2:
            graph_to_contig_map[n] = n2

    contig_graph = ContigGraph(
        graph=graph,
        path=gfa_file,
        contig_names=contig_names,
        contig_ids=None,
        contig_sequences=contig_sequences,
        contig_descriptions=contig_descriptions,
        graph_to_contig_map=graph_to_contig_map,
        self_loops=self_loops,
    )

    return contig_graph
