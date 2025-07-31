#!/usr/bin/env python3

import warnings

import pandas as pd


class ContigGraph:
    """
    Represents a contig-level assembly graph derived from a GFA file.

    This class encapsulates structural and sequence metadata for contigs constructed
    from GFA segment links, and optionally includes sequence, description, and
    graph-contig mappings.

    Attributes
    ----------
    graph : igraph.Graph
        The undirected graph representing the contig-level assembly graph.
    vcount : int
        The number of vertices in the graph.
    ecount : int
        The number of edges in the graph.
    file_path : str
        Path to the GFA file.
    contig_names : bidict
        Mapping from internal node IDs (starting from 0) to contig name.
    contig_parser : FastaParser
        FastaParser object containing the file pointers to contig sequences
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
        vcount,
        ecount,
        file_path,
        contig_names,
        contig_parser,
        contig_descriptions=None,
        graph_to_contig_map=None,
        self_loops=None,
    ):
        self.graph = graph
        self.vcount = vcount
        self.ecount = ecount
        self.file_path = file_path
        self.contig_names = contig_names  # node_id -> contig_name
        self.contig_parser = contig_parser
        self.contig_descriptions = (
            contig_descriptions  # name in contigs.fa (for MEGAHIT)
        )
        self.graph_to_contig_map = (
            graph_to_contig_map  # graph name -> contig.fa name (for MEGAHIT)
        )
        self.self_loops = self_loops

    def get_neighbours(self, contig_id: str) -> list:
        """
        Get neighbor contig IDs connected to the given contig.

        Parameters
        ----------
        contig_id : str
            The contig ID.

        Returns
        -------
        list of str
            List of neighboring contig IDs.
        """
        contig_names_rev = self.contig_names.inverse
        vid = contig_names_rev[contig_id]
        neighbor_ids = self.graph.neighbors(vid)
        return [self.contig_names[nid] for nid in neighbor_ids]

    def is_connected(self, from_contig: str, to_contig: str) -> bool:
        """
        Check if there is a path between two contigs in the graph.

        This method determines whether a path exists between the contig
        specified by `from_contig` and the contig specified by `to_contig`
        using the underlying graph's shortest path search.

        Parameters
        ----------
        from_contig : str
            Name of the starting contig.
        to_contig : str
            Name of the target contig.

        Returns
        -------
        bool
            True if there is a path connecting `from_contig` to `to_contig`,
            False otherwise.
        """
        contig_names_rev = self.contig_names.inverse
        from_id = contig_names_rev[from_contig]
        to_id = contig_names_rev[to_contig]

        with warnings.catch_warnings():
            # Suppress igraph's "RuntimeWarning: Couldn't reach some vertices"
            warnings.simplefilter("ignore")
            results = self.graph.get_shortest_paths(from_id, to=to_id)

        if len(results[0]) > 0:
            return True
        else:
            return False

    def get_adjacency_matrix(self, type="matrix"):
        """
        Return the adjacency matrix of the contig graph in different formats.

        Parameters
        ----------
        type : str, optional
            The return type. Options are:
            - "matrix": Return the adjacency matrix object from `self.graph.get_adjacency()`.
            - "pandas": Return a Pandas DataFrame with contig names as row and column labels.

        Returns
        -------
        adjacency : object or pandas.DataFrame
            - If `type="matrix"`, returns the adjacency matrix object.
            - If `type="pandas"`, returns a DataFrame where both rows and columns are indexed by contig names.

        Raises
        ------
        ValueError
            If `type` is not "matrix" or "pandas".
        """

        adj = self.graph.get_adjacency()

        if type == "matrix":
            return adj
        elif type == "pandas":
            labels = list(self.contig_names.values())
            adj_df = pd.DataFrame(adj, index=labels, columns=labels)
            return adj_df
        else:
            raise ValueError("type must be 'matrix' or 'pandas'")

    def get_connected_components(self) -> list:
        """
        Calculate the average node degree of the graph.

        Returns
        -------
        list
            A list of the connected components
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
        """
        contig_lengths = [
            len(self.contig_parser.get_sequence(seq))
            for seq in self.contig_names.values()
        ]
        return sum(contig_lengths)

    def calculate_average_contig_length(self) -> int:
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
        """

        contig_lengths = [
            len(self.contig_parser.get_sequence(seq))
            for seq in self.contig_names.values()
        ]
        if len(contig_lengths) == 0:
            raise ValueError(
                "Graph does not have any segments, cannot calculate average segment length"
            )

        return int(sum(contig_lengths) / len(contig_lengths))

    def calculate_n50_l50(self) -> tuple[int, int]:
        """
        Calculate N50 and L50 from a list of segment lengths.

        Returns
        -------
        tuple of (int, int)
            A tuple containing:
            - N50 : int
                The length N such that 50% of the total length is contained in segments of length ≥ N.
            - L50 : int
                The minimum number of segments whose summed length ≥ 50% of the total.
        """

        contig_lengths = [
            len(self.contig_parser.get_sequence(seq))
            for seq in self.contig_names.values()
        ]
        sorted_lengths = sorted(contig_lengths, reverse=True)
        total_length = sum(sorted_lengths)
        cum_sum = 0

        for i, length in enumerate(sorted_lengths):
            cum_sum += length
            if cum_sum >= total_length / 2:
                return length, i + 1

    def get_gc_content(self) -> float:
        """
        Calculate the GC content of sequences.

        Returns
        -------
        float
            GC content as a percentage of total base pairs.

        Raises
        ------
        ValueError
            If total length of the segments is zero.
        """

        contig_sequences = [
            self.contig_parser.get_sequence(seq) for seq in self.contig_names.values()
        ]
        total_length = self.calculate_total_length()

        if total_length == 0:
            raise ValueError(
                "Total length of segments is zero, cannot calculate GC content"
            )

        gc_count = sum(seq.count("G") + seq.count("C") for seq in contig_sequences)
        return gc_count / total_length
