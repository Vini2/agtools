#!/usr/bin/env python3

import io

from bidict import bidict
from Bio import SeqIO
from igraph import Graph

from agtools.core.contig_graph import ContigGraph
from agtools.core.fasta_parser import FastaParser


def _get_links_and_contig_mapping_myloasm(gfa_file: str) -> tuple:
    """
    Parse a GFA file to extract contig information and connectivity
    information (links) between contigs.

    Parameters
    ----------
    gfa_file : str
        Path to the MEGAHIT-style GFA file.

    Returns
    -------
    tuple
    contig_names : list
        List of contig names
    contig_name_to_id : dict
        Mapping from contig name to internal ID
    graph_contig_seqs : dict
        Mapping of segment ID -> sequence length in graph file.
    edge_list : list
        List of edges
    self_loops : list
        List of self loops
    lcount : int
        The number of links (lines starting with tag "L") in the graph
    """
    lcount = 0
    contig_names = list()
    contig_name_to_id = dict()
    edge_list = set()
    self_loops = set()

    graph_contig_seqs = {}

    with io.open(gfa_file, mode="r", buffering=1024 * 1024) as f:
        while True:
            line = f.readline()
            if not line:
                break

            tag = line[0]

            if not line:
                continue

            if tag == "S":  # Segment line
                parts = line.rstrip().split("\t")
                contig_name = parts[1]
                seq = parts[2]

                contig_id = len(contig_names)
                contig_name_to_id[contig_name] = contig_id
                contig_names.append(contig_name)

                graph_contig_seqs[contig_name] = len(seq)

    with io.open(gfa_file, mode="r", buffering=1024 * 1024) as f:
        while True:
            line = f.readline()
            if not line:
                break

            tag = line[0]

            if not line:
                continue

            if tag == "L":  # Link line
                lcount += 1
                parts = line.rstrip().split("\t")
                from_seg, from_orient = parts[1], parts[2]
                to_seg, to_orient = parts[3], parts[4]
                overlap = int(parts[5][:-1])  # Remove trailing M

                source = contig_name_to_id[from_seg]
                target = contig_name_to_id[to_seg]

                if source == target:
                    self_loops.add(source)
                else:
                    edge_list.add((source, target))

    return (
        contig_names,
        contig_name_to_id,
        graph_contig_seqs,
        list(edge_list),
        list(self_loops),
        lcount,
    )


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

    (
        contig_names,
        contig_name_to_id,
        graph_contig_seqs,
        edge_list,
        self_loops,
        lcount,
    ) = _get_links_and_contig_mapping_myloasm(gfa_file)

    # Create graph
    graph = Graph()

    # Add vertices
    graph.add_vertices(len(contig_names))

    # Name vertices with contig identifiers
    graph.vs["label"] = contig_names

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
        lcount=lcount,
        ecount=graph.ecount(),
        file_path=gfa_file,
        contig_names=contig_names,
        contig_name_to_id=contig_name_to_id,
        contig_parser=parser,
        contig_descriptions=contig_descriptions,
        graph_to_contig_map=graph_to_contig_map,
        self_loops=self_loops,
    )

    return contig_graph
