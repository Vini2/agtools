#!/usr/bin/env python3

"""
Graph utility classes and functions for parsing and manipulating assembly graphs.

Includes:
- `UnitigGraph`: Represents unitig-level assembly graphs from GFA format.
- `ContigGraph`: Represents contig-level graphs with optional mappings.
- `parse_fastg`: Parses FASTG files into segments and adjacency data.

Dependencies:
- `igraph` for graph structures
- `bidict` for bidirectional maps
- `Bio.Seq.Seq` for sequence storage
"""

import re
from collections import defaultdict

from bidict import bidict
from Bio.Seq import Seq
from igraph import Graph


class UnitigGraph:
    """
    Represents an assembly graph parsed from a GFA file.

    Attributes
    ----------
    graph : igraph.Graph
        The undirected graph representing the unitig-level assembly graph.
    path : str
        Path to the original GFA file.
    oriented_links : dict
        Mapping from (from_seg, to_seg) → list of (from_orient, to_orient).
    link_overlap : dict
        Mapping from oriented segment pair to overlap length.
    segment_names : bidict
        Maps internal node IDs to Segment IDs.
    segment_names_rev : bidict
        Reverse mapping of segment names to node IDs.
    segment_sequences : dict
        Segment ID → sequence (as Bio.Seq.Seq).
    segment_lengths : dict
        Segment ID → length of sequence.
    self_loops : list
        List of segment IDs that form self-loops.
    """

    def __init__(self):
        self.graph = Graph(directed=False)
        self.path = None
        self.oriented_links = defaultdict(lambda: defaultdict(list))
        self.link_overlap = dict()
        self.segment_names = bidict()  # node_id → segment_id
        self.segment_names_rev = None  # segment_id → node_id
        self.segment_sequences = dict()  # segment_id → sequence
        self.segment_lengths = dict()  # segment_id → length
        self.self_loops = []

    @classmethod
    def from_gfa(cls, path: str) -> "UnitigGraph":
        """
        Parse a GFA file into a UnitigGraph object.

        Parameters
        ----------
        path : str
            Path to the GFA file.

        Returns
        -------
        UnitigGraph
            The constructed graph object with segments, links, and metadata.
        """

        ug = cls()
        node_count = 0
        links = []

        ug.path = path

        with open(path) as f:
            for line in f:

                if line.startswith("S"):
                    parts = line.strip().split("\t")
                    seg_id = parts[1]
                    seq = parts[2]
                    ug.segment_names[node_count] = seg_id
                    ug.segment_sequences[seg_id] = Seq(seq)
                    ug.segment_lengths[seg_id] = len(seq)
                    node_count += 1
                elif line.startswith("L"):
                    parts = line.strip().split("\t")
                    from_seg, from_orient = parts[1], parts[2]
                    to_seg, to_orient = parts[3], parts[4]
                    overlap = int(parts[5][:-1])  # Remove trailing M

                    links.append((from_seg, to_seg))
                    ug._add_oriented_links(
                        from_seg, to_seg, from_orient, to_orient, overlap
                    )

        ug.segment_names_rev = ug.segment_names.inverse
        ug.graph.add_vertices(node_count)

        for i in range(node_count):
            seg_name = ug.segment_names[i]
            ug.graph.vs[i]["id"] = i
            ug.graph.vs[i]["name"] = seg_name
            ug.graph.vs[i]["label"] = f"{seg_name}\nID:{i}"

        edge_list, ug.self_loops = ug._get_graph_edges(links)
        ug.graph.add_edges(edge_list)
        ug.graph.simplify(multiple=True, loops=False, combine_edges=None)
        return ug

    def _add_oriented_links(self, from_seg, to_seg, from_orient, to_orient, overlap):
        """
        Store oriented link and its overlap, along with its symmetric reverse.

        Parameters
        ----------
        from_seg : str
            Source segment ID.
        to_seg : str
            Destination segment ID.
        from_orient : str
            Orientation of source ('+' or '-').
        to_orient : str
            Orientation of destination.
        overlap : int
            Overlap length in base pairs.
        """

        key1 = f"{from_seg}{from_orient}"
        key2 = f"{to_seg}{to_orient}"
        self.oriented_links[from_seg][to_seg].append((from_orient, to_orient))
        self.link_overlap[(key1, key2)] = overlap

        # Add symmetric reverse
        rev1 = "+" if from_orient == "-" else "-"
        rev2 = "+" if to_orient == "-" else "-"
        self.oriented_links[to_seg][from_seg].append((rev2, rev1))
        self.link_overlap[(f"{to_seg}{rev2}", f"{from_seg}{rev1}")] = overlap

    def _get_graph_edges(self, links):
        """
        Convert parsed segment links into edge list for igraph.

        Parameters
        ----------
        links : list of tuple
            Pairs of segment IDs representing links.

        Returns
        -------
        tuple
            (edges, self_loops) where edges is a list of (src_id, tgt_id),
            and self_loops is a list of segment IDs forming loops.
        """

        edges = []
        loops = []
        for from_edge, to_edge in links:
            if from_edge == to_edge:
                loops.append(from_edge)
            else:
                src = self.segment_names_rev[from_edge]
                tgt = self.segment_names_rev[to_edge]
                edges.append((src, tgt))
        return edges, loops

    def get_neighbors(self, seg_id: str) -> list:
        """
        Get neighbor segment IDs connected to the given segment.

        Parameters
        ----------
        seg_id : str
            The segment ID.

        Returns
        -------
        list of str
            List of neighboring segment IDs.
        """

        vid = self.segment_names_rev[seg_id]
        neighbor_ids = self.graph.neighbors(vid)
        return [self.segment_names[nid] for nid in neighbor_ids]


class ContigGraph:
    """
    Represents a contig-level graph derived from an assembly graph.

    This class encapsulates structural and sequence metadata for contigs constructed
    from GFA segment links, and optionally includes sequence, description, and
    graph-contig mappings.

    Attributes
    ----------
    graph : igraph.Graph
        The igraph object representing the contig-level graph structure.
    path : str
        Path to the GFA file.
    contig_names : bidict
        Mapping from node ID to contig name string.
    contig_ids : bidict, optional
        Mapping from node ID to internal or external contig number.
    contig_sequences : dict[str, str], optional
        Dictionary mapping contig names to DNA sequences.
    contig_descriptions : dict[str, str], optional
        Dictionary mapping contig names to additional descriptions in FASTA file.
    graph_to_contig_map : dict[int, str], optional
        Dictionary mapping from unitig-level node IDs to contig identifiers
    self_loops : list[str], optional
        List of contig names that form self-loops in the graph.
    """

    def __init__(
        self,
        graph,
        path,
        contig_names,
        contig_ids=None,
        contig_sequences=None,
        contig_descriptions=None,
        graph_to_contig_map=None,
        self_loops=None,
    ):
        self.graph = graph
        self.path = path
        self.contig_names = contig_names  # node_id → segment_id
        self.contig_ids = contig_ids  # node_id → contig_i
        self.contig_sequences = contig_sequences
        self.contig_descriptions = contig_descriptions
        self.graph_to_contig_map = graph_to_contig_map  # for MEGAHIT
        self.self_loops = self_loops


def parse_fastg(fastg_file: str) -> tuple:
    """
    Parse a FASTG file and extract segment sequences and edges.

    Parameters
    ----------
    fastg_file : str
        Path to the FASTG file.

    Returns
    -------
    tuple
        segments : dict
            Mapping from node ID to DNA sequence.
        edges : dict
            Mapping from node ID to list of adjacent node IDs.
    """

    segments = {}
    edges = {}

    with open(fastg_file, "r") as f:
        current_node = None
        sequence = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_node and sequence:
                    segments[current_node] = "".join(sequence)
                    sequence = []

                header = line[1:]
                parts = header.split(":")
                node = parts[0].strip("'")
                neighbors = []
                if len(parts) > 1:
                    neighbors = re.split(r"[,\s]+", parts[1])
                current_node = node
                edges[node] = neighbors
            else:
                sequence.append(line)

        if current_node and sequence:
            segments[current_node] = "".join(sequence)

    return segments, edges
