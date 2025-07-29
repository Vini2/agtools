#!/usr/bin/env python3

from collections import defaultdict

import pandas as pd
from bidict import bidict
from Bio.Seq import Seq
from igraph import Graph


class UnitigGraph:
    """
    Represents a unitig-level assembly graph parsed from a GFA file.

    Attributes
    ----------
    graph : igraph.Graph
        The undirected graph representing the unitig-level assembly graph.
    vcount : int
        The number of vertices in the graph.
    ecount : int
        The number of edges in the graph.
    file_path : str
        Path to the original GFA file.
    oriented_links : dict
        Mapping from (from_seg, to_seg) -> list of (from_orient, to_orient).
    link_overlap : dict
        Mapping from oriented segment pair to overlap length.
    segment_names : bidict
        Maps internal node IDs (starting from 0) to Segment IDs.
    segment_lengths : dict
        Segment ID -> length of sequence.
    segment_offsets : dict
        Segment ID -> byte offset to the segment line in the gfa file.
    self_loops : list
        List of segment IDs that form self-loops.
    self.paths :

    References
    ----------
    GFA: Graphical Fragment Assembly (GFA) Format Specification
    https://github.com/GFA-spec/GFA-spec
    """

    def __init__(self):
        self.graph = Graph(directed=False)
        self.vcount = 0
        self.ecount = 0
        self.file_path = None
        self.oriented_links = defaultdict(lambda: defaultdict(list))
        self.link_overlap = dict()
        self.segment_names = bidict()  # node_id -> segment_name
        self.segment_sequences = dict()  # segment_id -> sequence
        self.segment_lengths = dict()  # segment_id -> length
        self.segment_offsets = dict()  # segment_id -> byte offset in file
        self.self_loops = []

    @classmethod
    def from_gfa(cls, file_path: str) -> "UnitigGraph":
        """
        Parse a GFA file into a UnitigGraph object.

        Parameters
        ----------
        file_path : str
            Path to the GFA file.

        Returns
        -------
        UnitigGraph
            The constructed graph object with segments, links, and metadata.
        """

        ug = cls()
        node_count = 0
        links = []

        ug.file_path = file_path

        with open(file_path) as f:
            while True:
                pos = f.tell()
                line = f.readline()
                if not line:
                    break

                if line.startswith("S"):
                    parts = line.strip().split("\t")
                    seg_name = parts[1]
                    seq = parts[2]
                    ug.segment_names[node_count] = seg_name
                    ug.segment_offsets[seg_name] = pos
                    ug.segment_lengths[seg_name] = len(seq)
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

        ug.graph.add_vertices(node_count)

        for i in range(node_count):
            seg_name = ug.segment_names[i]
            ug.graph.vs[i]["id"] = i
            ug.graph.vs[i]["name"] = seg_name

        edge_list, ug.self_loops = ug._get_graph_edges(links)
        ug.graph.add_edges(edge_list)
        ug.graph.simplify(multiple=True, loops=False, combine_edges=None)

        ug.vcount = ug.graph.vcount()
        ug.ecount = ug.graph.ecount()

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
        segment_names_rev = self.segment_names.inverse
        for from_edge, to_edge in links:
            if from_edge == to_edge:
                loops.append(from_edge)
            else:
                src = segment_names_rev[from_edge]
                tgt = segment_names_rev[to_edge]
                edges.append((src, tgt))
        return edges, loops

    def get_segment_sequence(self, seg_name: str) -> Seq:
        """
        Load the DNA sequence for a specific segment on demand.

        This method retrieves the sequence of a segment from the original GFA file
        using byte offsets, without loading all sequences into memory at once.

        Parameters
        ----------
        seg_name : str
            The segment identifier (ID) whose DNA sequence should be retrieved.

        Returns
        -------
        Bio.Seq.Seq
            The DNA sequence corresponding to the given segment.

        Raises
        ------
        KeyError
            If the segment name does not exist in the graph.
        ValueError
            If the retrieved sequence length does not match the expected length
            recorded during graph construction.
        """
        pos = self.segment_offsets[seg_name]
        with open(self.file_path, "r") as f:
            f.seek(pos)
            line = f.readline()
            seq = line.strip().split("\t")[2]

            if seg_name not in self.segment_lengths:
                raise KeyError("Segment name does not exist in the graph")
            if len(seq) == self.segment_lengths[seg_name]:
                return Seq(seq)
            else:
                raise ValueError("Wrong sequence retrieved")

    def get_neighbours(self, seg_id: str) -> list:
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
        segment_names_rev = self.segment_names.inverse
        vid = segment_names_rev[seg_id]
        neighbor_ids = self.graph.neighbors(vid)
        return [self.segment_names[nid] for nid in neighbor_ids]

    def get_adjacency_matrix(self, type="matrix"):
        """
        Return the adjacency matrix of the unitig graph in different formats.

        Parameters
        ----------
        type : str, optional
            The return type. Options are:
            - "matrix": Return the adjacency matrix object from `self.graph.get_adjacency()`.
            - "pandas": Return a Pandas DataFrame with unitig names as row and column labels.

        Returns
        -------
        adjacency : object or pandas.DataFrame
            - If `type="matrix"`, returns the adjacency matrix object.
            - If `type="pandas"`, returns a DataFrame where both rows and columns are indexed by unitig names.

        Raises
        ------
        ValueError
            If `type` is not "matrix" or "pandas".
        """

        adj = self.graph.get_adjacency()

        if type == "matrix":
            return adj
        elif type == "pandas":
            labels = list(self.segment_names.values())
            adj_df = pd.DataFrame(adj, index=labels, columns=labels)
            return adj_df
        else:
            raise ValueError("type must be 'matrix' or 'pandas'")

    def is_connected(self, from_seg: str, to_seg: str) -> bool:
        """
        Check if there is a path between two segments in the graph.

        This method determines whether a path exists between the segment
        specified by `from_seg` and the segment specified by `to_seg`
        using the underlying graph's shortest path search.

        Parameters
        ----------
        from_seg : str
            Name of the starting segment.
        to_seg : str
            Name of the target segment.

        Returns
        -------
        bool
            True if there is a path connecting `from_seg` to `to_seg`,
            False otherwise.
        """
        segments_names_rev = self.segment_names.inverse
        from_id = segments_names_rev[from_seg]
        to_id = segments_names_rev[to_seg]

        results = self.graph.get_shortest_paths(from_id, to=to_id)

        if len(results[0]) > 0:
            return True
        else:
            return False
