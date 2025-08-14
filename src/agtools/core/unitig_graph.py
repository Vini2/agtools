#!/usr/bin/env python3

from collections import defaultdict

import pandas as pd
from Bio.Seq import Seq
from igraph import Graph


class UnitigGraph:
    """
    Represents a unitig-level assembly graph parsed from a GFA file.

    Attributes
    ----------
    graph : igraph.Graph
        The undirected graph representing the unitig-level assembly graph
    vcount : int
        The number of vertices (segments) in the graph
    lcount : int
        The number of links (lines starting with tag "L") in the graph
    ecount : int
        The number of edges in the graph after simplification
    pcount : int
        The number of paths (lines starting with tag "P") in the graph
    file_path : str
        Path to the original GFA file
    segment_names : list
        List of segment names
    segment_name_to_id : dict
        Mapping from segment name to internal ID
        This is used to map segment names to their vertex IDs in the graph
    segment_lengths : dict
        Mapping from segment name to length of sequence
    segment_offsets : dict
        Mapping from segment name to byte offset of the segment line in the gfa file
    oriented_links : dict
        Mapping from [from segment id][to segment id] -> list of (from orientation, to orientation)
    link_overlap : dict
        Mapping from oriented segment pair (from segment id, from orientation, to segment id, to orientarion) -> overlap length
    paths: dict
        Mapping from path name -> (segment names, overlaps)
    self_loops : list
        List of segment IDs that form self-loops.

    Methods
    -------
    from_gfa(file_path)
        Parse a GFA file into a UnitigGraph object.
    get_segment_sequence(seg_name)
        Retrieve a DNA sequence for a segment.
    get_neighbors(seg_id)
        Get neighboring segments of a given segment.
    get_adjacency_matrix(type="matrix")
        Return the adjacency matrix as igraph or pandas DataFrame.
    is_connected(from_seg, to_seg)
        Check if there is a path between two segments in the graph.
    get_connected_components()
        Get connected components of the graph.
    calculate_average_node_degree()
        Calculate the average node degree of the graph.
    calculate_total_length()
        Calculate the total length of all segments in the graph.
    calculate_average_segment_length()
        Calculate the average segment length.
    calculate_n50_l50()
        Calculate N50 and L50 for the segments in the graph.
    get_gc_content()
        Calculate the GC content of segment sequences.

    Examples
    --------
    >>> from agtools.core.unitig_graph import UnitigGraph
    >>> ug = UnitigGraph.from_gfa("assembly.gfa")
    >>> ug.vcount
    42
    >>> ug.ecount
    80

    References
    ----------
    GFA: Graphical Fragment Assembly (GFA) Format Specification
    https://github.com/GFA-spec/GFA-spec
    """

    __slots__ = (
        "graph",
        "vcount",
        "lcount",
        "ecount",
        "pcount",
        "file_path",
        "oriented_links",
        "link_overlap",
        "paths",
        "segment_names",
        "segment_name_to_id",
        "segment_lengths",
        "segment_offsets",
        "self_loops",
    )

    def __init__(self):
        self.graph = Graph(directed=False)
        self.vcount = 0
        self.lcount = 0
        self.ecount = 0
        self.pcount = 0
        self.file_path = None
        self.segment_names = list()         # list of segment names
        self.segment_name_to_id = dict()    # segment name -> internal ID
        self.segment_lengths = dict()       # segment_id -> length
        self.segment_offsets = dict()       # segment_id -> byte offset in file
        self.oriented_links = defaultdict(lambda: defaultdict(list))
        self.link_overlap = dict()
        self.paths = dict()                 # path_id -> segment names
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

        Examples
        --------
        >>> ug = UnitigGraph.from_gfa("assembly.gfa")
        >>> ug.vcount
        42
        >>> ug.ecount
        80
        """

        ug = cls()
        edge_list = []

        ug.file_path = file_path

        with open(file_path) as f:
            while True:
                pos = f.tell()
                line = f.readline()
                if not line:
                    break

                tag = line[0]

                if not line:
                    continue
                
                if tag == "S":      # Segment line
                    parts = line.rstrip().split("\t")
                    seg_name = parts[1]
                    seq = parts[2]
                    seg_id = len(ug.segment_names)
                    ug.segment_name_to_id[seg_name] = seg_id
                    ug.segment_names.append(seg_name)
                    ug.segment_offsets[seg_name] = pos
                    ug.segment_lengths[seg_name] = len(seq)

                elif tag == "L":    # Link line
                    ug.lcount += 1
                    parts = line.rstrip().split("\t")
                    from_seg, from_orient = parts[1], parts[2]
                    to_seg, to_orient = parts[3], parts[4]
                    overlap = int(parts[5][:-1])  # Remove trailing M

                    source = ug.segment_name_to_id[from_seg]
                    target = ug.segment_name_to_id[to_seg]

                    if source == target:
                        ug.self_loops.append(source)
                    else:
                        edge_list.append((source, target))

                    ug.oriented_links[source][target].append((from_orient, to_orient))
                    ug.link_overlap[(source, from_orient, target, to_orient)] = overlap

                    # Add symmetric reverse
                    rev1 = "+" if from_orient == "-" else "-"
                    rev2 = "+" if to_orient == "-" else "-"
                    ug.oriented_links[target][source].append((rev2, rev1))
                    ug.link_overlap[(target, rev2, source, rev1)] = overlap
                    
                elif tag == "P":    # Path line
                    ug.pcount += 1
                    parts = line.rstrip().split("\t")
                    path_name = parts[1]
                    segment_tokens = parts[2].split(",")
                    # segment_ids = [f"{ug.segment_name_to_id[segment[:-1]]}{segment[-1]}" for segment in segment_tokens]
                    overlaps = [int(x[:-1]) if x.endswith(("M", "m")) else x for x in parts[3].split(",")]
                    ug.paths[path_name] = (segment_tokens, overlaps)

        # Add vertices
        ug.vcount = len(ug.segment_names)
        ug.graph.add_vertices(ug.vcount)
        ug.graph.vs["name"] = ug.segment_names

        # Add edges
        ug.graph.add_edges(edge_list)
        ug.graph.simplify(multiple=True, loops=False, combine_edges=None)
        ug.ecount = ug.graph.ecount()

        return ug
    

    def get_segment_sequence(self, seg_name: str) -> Seq:
        """
        Retrieve a DNA sequence for a segment.

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

        Examples
        --------
        >>> ug.get_segment_sequence("unitig_1")[:10]
        Seq('ATGCGTACGG')
        """
        pos = self.segment_offsets[seg_name]
        with open(self.file_path, "r") as f:
            f.seek(pos)
            line = f.readline()
            seq = line.rstrip().split("\t")[2]

            if seg_name not in self.segment_lengths:
                raise KeyError("Segment name does not exist in the graph")
            if len(seq) == self.segment_lengths[seg_name]:
                return Seq(seq)
            else:
                raise ValueError("Wrong sequence retrieved")

    def get_neighbors(self, seg_id: str) -> list:
        """
        Get neighboring segments of a given segment.

        Parameters
        ----------
        seg_id : str
            The segment ID.

        Returns
        -------
        list of str
            List of neighboring segment IDs.

        Examples
        --------
        >>> ug.get_neighbors("unitig_1")
        ['unitig_2', 'unitig_3']
        """
        vid = self.segment_name_to_id[seg_id]
        return [self.segment_names[nid] for nid in self.graph.neighbors(vid)]

    def get_adjacency_matrix(self, type="matrix"):
        """
        Return the adjacency matrix as igraph or pandas DataFrame.

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

        Examples
        --------
        >>> matrix = ug.get_adjacency_matrix()
        >>> isinstance(matrix, list)
        True
        >>> df = ug.get_adjacency_matrix(type="pandas")
        >>> df.head()
                    unitig_1  unitig_2  unitig_3
        unitig_1          0         1         0
        unitig_2          1         0         1
        unitig_3          0         1         0
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

        Examples
        --------
        >>> ug.is_connected("unitig_1", "unitig_2")
        True
        """
        from_id = self.segment_name_to_id[from_seg]
        to_id = self.segment_name_to_id[to_seg]

        results = self.graph.get_shortest_paths(from_id, to=to_id)

        if len(results[0]) > 0:
            return True
        else:
            return False

    def get_connected_components(self) -> list:
        """
        Get connected components of the graph.

        Returns
        -------
        list
            A list of the connected components

        Examples
        --------
        >>> components = ug.get_connected_components()
        >>> len(components)
        3
        >>> [len(c) for c in components]
        [10, 8, 5]
        >>> components[0]
        [0, 1, 2, 3, ...]
        """
        return self.graph.components()

    def calculate_average_node_degree(self) -> int:
        """
        Calculate the average node degree of the graph.

        Returns
        -------
        int
            Average node degree of the graph.

        Raises
        ------
        ValueError
            If the graph does not have any segments.

        Examples
        --------
        >>> ug.calculate_average_node_degree()
        2
        """

        if self.graph.vcount() == 0:
            raise ValueError(
                "Graph does not have any segments, cannot calculate average node degree"
            )

        return int(sum(self.graph.degree()) / self.graph.vcount())

    def calculate_total_length(self) -> int:
        """
        Calculate the total length of all segments in the graph.

        Returns
        -------
        int
            Total length of all segments.

        Examples
        --------
        >>> ug.calculate_total_length()
        350000
        """
        return sum(self.segment_lengths.values())

    def calculate_average_segment_length(self) -> int:
        """
        Calculate the average segment length.

        Returns
        -------
        int
            Average segment length.

        Raises
        ------
        ValueError
            If the graph does not have any segments.

        Examples
        --------
        >>> ug.calculate_average_segment_length()
        3494
        """

        segment_lengths = self.segment_lengths
        if len(segment_lengths) == 0:
            raise ValueError(
                "Graph does not have any segments, cannot calculate average segment length"
            )

        return int(sum(segment_lengths.values()) / len(segment_lengths))

    def calculate_n50_l50(self) -> tuple[int, int]:
        """
        Calculate N50 and L50 for the segment in the graph.

        Returns
        -------
        tuple of (int, int)
            A tuple containing:
            - N50 : int
                The length N such that 50% of the total length is contained in segments of length ≥ N.
            - L50 : int
                The minimum number of segments whose summed length ≥ 50% of the total.

        Examples
        --------
        >>> ug.calculate_n50_l50()
        (15000, 12)
        """

        lengths = self.segment_lengths.values()
        sorted_lengths = sorted(lengths, reverse=True)
        total_length = sum(sorted_lengths)
        cum_sum = 0

        for i, length in enumerate(sorted_lengths):
            cum_sum += length
            if cum_sum >= total_length / 2:
                return length, i + 1

    def get_gc_content(self) -> float:
        """
        Calculate the GC content of segment sequences.

        Returns
        -------
        float
            GC content as a percentage of total base pairs.

        Raises
        ------
        ValueError
            If total length of the segments is zero.

        Examples
        --------
        >>> round(ug.get_gc_content(), 2)
        0.42
        """

        sequences = [
            self.get_segment_sequence(seq) for seq in self.segment_lengths.keys()
        ]
        total_length = self.calculate_total_length()

        if total_length == 0:
            raise ValueError(
                "Total length of segments is zero, cannot calculate GC content"
            )

        gc_count = sum(seq.count("G") + seq.count("C") for seq in sequences)
        return gc_count / total_length
