#!/usr/bin/env python3

import io

from igraph import Graph

from agtools.core.contig_graph import ContigGraph
from agtools.core.fasta_parser import FastaParser


def _get_links_and_contig_mapping_myloasm(gfa_file: str, contig_index: dict) -> tuple:
    """
    Parse a GFA file to extract contig information and connectivity
    information (links) between contigs.

    Parameters
    ----------
    gfa_file : str
        Path to the myloasm-style GFA file.
    contig_index: dict
        Dictionary of contigs in the FASTA file.

    Returns
    -------
    tuple
    contig_names : list
        List of contig names
    contig_name_to_id : dict
        Mapping from contig name to internal ID
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

                if contig_name in contig_index:
                    contig_id = len(contig_names)
                    contig_name_to_id[contig_name] = contig_id
                    contig_names.append(contig_name)

            elif tag == "L":  # Link line
                lcount += 1
                parts = line.rstrip().split("\t")
                from_seg, from_orient = parts[1], parts[2]
                to_seg, to_orient = parts[3], parts[4]
                overlap = int(parts[5][:-1])  # Remove trailing M

                if from_seg in contig_index and to_seg in contig_index:

                    source = contig_name_to_id[from_seg]
                    target = contig_name_to_id[to_seg]

                    if source == target:
                        self_loops.add(source)
                    else:
                        edge_list.add((source, target))

    return contig_names, contig_name_to_id, list(edge_list), list(self_loops), lcount


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
    contig_names, contig_name_to_id, edge_list, self_loops, lcount = (
        _get_links_and_contig_mapping_myloasm(gfa_file, parser.index)
    )

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

    # Create ContigGraph object
    contig_graph = ContigGraph(
        graph=graph,
        vcount=graph.vcount(),
        lcount=lcount,
        ecount=graph.ecount(),
        file_path=gfa_file,
        contig_names=contig_names,
        contig_name_to_id=contig_name_to_id,
        contig_parser=parser,
        contig_descriptions=None,
        graph_to_contig_map=None,
        self_loops=self_loops,
    )

    return contig_graph
