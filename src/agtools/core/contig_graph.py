#!/usr/bin/env python3

import warnings

import pandas as pd
from Bio.Seq import Seq


class ContigGraph:
    """
    Represents a contig-level assembly graph derived from a GFA file.

    Attributes
    ----------
    graph : igraph.Graph
        The undirected graph representing the contig-level assembly graph.
    vcount : int
        The number of vertices (contigs) in the graph.
    lcount : int
        The number of links (lines starting with tag "L") in the GFA file.
    ecount : int
        The number of edges in the graph after simplification
    file_path : str
        Path to the GFA file.
    contig_names : list
        List of contig names.
    contig_name_to_id : dict
        Mapping from contig name to internal ID.
        This is used to map contig names to their vertex IDs in the graph.
    contig_parser : FastaParser
        FastaParser object containing the file pointers to contig sequences.
    contig_descriptions : dict[str, str], optional
        Dictionary mapping contig names to additional descriptions in FASTA file.
    graph_to_contig_map : bidict[int, str], optional
        Bi-directional dictionary mapping from contig identifiers in the GFA file to FASTA file.
    self_loops : list[str], optional
        List of contig names that form self-loops in the graph.

    Methods
    -------
    get_contig_sequence(contig_id)
        Retrieve a DNA sequence for a contig.
    get_neighbors(contig_name)
        Get neighboring contigs of a given contig.
    get_adjacency_matrix(type="matrix")
        Return the adjacency matrix as igraph or pandas DataFrame.
    is_connected(from_contig, to_contig)
        Check if there is a path between two contigs in the graph.
    get_connected_components()
        Get connected components of the graph.
    calculate_average_node_degree()
        Calculate the average node degree of the graph.
    calculate_total_length()
        Calculate the total length of all contigs in the graph.
    calculate_average_contig_length()
        Calculate the average contig length.
    calculate_n50_l50()
        Calculate N50 and L50 for the contigs in the graph.
    get_gc_content()
        Calculate the GC content of contig sequences.
    """

    __slots__ = (
        "graph",
        "vcount",
        "lcount",
        "ecount",
        "file_path",
        "contig_names",
        "contig_name_to_id",
        "contig_parser",
        "contig_descriptions",
        "graph_to_contig_map",
        "self_loops",
    )

    def __init__(
        self,
        graph,
        vcount,
        lcount,
        ecount,
        file_path,
        contig_names,
        contig_name_to_id,
        contig_parser,
        contig_descriptions=None,
        graph_to_contig_map=None,
        self_loops=None,
    ):
        self.graph = graph
        self.vcount = vcount
        self.lcount = lcount
        self.ecount = ecount
        self.file_path = file_path
        self.contig_names = contig_names  # list of contig names
        self.contig_name_to_id = contig_name_to_id  # contig_name -> node_id
        self.contig_parser = contig_parser
        self.contig_descriptions = (
            contig_descriptions  # name in contigs.fa (for MEGAHIT)
        )
        self.graph_to_contig_map = (
            graph_to_contig_map  # graph name -> contig.fa name (for MEGAHIT)
        )
        self.self_loops = self_loops

    def get_contig_sequence(self, contig_name: str) -> Seq:
        """
        Retrieve a DNA sequence for a contig.

        This method retrieves the sequence of a contig from the contigs file
        using byte offsets, without loading all sequences into memory at once.

        Parameters
        ----------
        contig_name : str
            The contig identifier whose DNA sequence should be retrieved.

        Returns
        -------
        Bio.Seq.Seq
            The DNA sequence corresponding to the given contig.

        Examples
        --------
        >>> cg.get_contig_sequence("contig_1")
        Seq('TTGATGCGACGTACGG')
        """
        return self.contig_parser.get_sequence(contig_name)

    def get_neighbors(self, contig_name: str) -> list:
        """
        Get neighboring contigs of a given contig.

        Parameters
        ----------
        contig_name : str
            The contig name.

        Returns
        -------
        list of str
            List of neighboring contig names.

        Examples
        --------
        >>> cg.get_neighbors("contig_1")
        ['contig_2', 'contig_3']
        """
        vid = self.contig_name_to_id[contig_name]
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

        Raises
        ------
        KeyError
            If the contig names do not exist in the assembly.

        Examples
        --------
        >>> cg.is_connected("contig_1", "contig_2")
        True
        """

        if (
            from_contig in self.contig_name_to_id
            and to_contig in self.contig_name_to_id
        ):
            from_id = self.contig_name_to_id[from_contig]
            to_id = self.contig_name_to_id[to_contig]

            with warnings.catch_warnings():
                # Suppress igraph's "RuntimeWarning: Couldn't reach some vertices"
                warnings.simplefilter("ignore")
                results = self.graph.get_shortest_paths(from_id, to=to_id)

            if len(results[0]) > 0:
                return True
            else:
                return False

        else:
            raise KeyError("Contig names do not exist in the assembly")

    def get_adjacency_matrix(self, type="matrix"):
        """
        Return the adjacency matrix as igraph or pandas DataFrame.

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

        Examples
        --------
        >>> matrix = cg.get_adjacency_matrix()
        >>> isinstance(matrix, list)
        True
        >>> df = cg.get_adjacency_matrix(type="pandas")
        >>> df.head()
                    contig_1  contig_2  contig_3
        contig_1          0         1         0
        contig_2          1         0         1
        contig_3          0         1         0
        """

        adj = self.graph.get_adjacency()

        if type == "matrix":
            return adj
        elif type == "pandas":
            labels = self.contig_names
            adj_df = pd.DataFrame(adj, index=labels, columns=labels)
            return adj_df
        else:
            raise ValueError("type must be 'matrix' or 'pandas'")

    def get_connected_components(self) -> list:
        """
        Get connected components of the graph.

        Returns
        -------
        list
            A list of the connected components with internal contig IDs.

        Examples
        --------
        >>> components = cg.get_connected_components()
        >>> len(components)
        3
        >>> [len(c) for c in components]
        [10, 8, 5]
        >>> components[0]
        [0, 1, 2, 3, ...]
        """
        return self.graph.components()

    def calculate_average_node_degree(self) -> float:
        """
        Calculate the average node degree of the graph.

        Returns
        -------
        int
            Average node degree of the graph.

        Raises
        ------
        ValueError
            If the graph does not have any contigs.

        Examples
        --------
        >>> cg.calculate_average_node_degree()
        1
        """

        if self.graph.vcount() == 0:
            raise ValueError(
                "Graph does not have any contigs, cannot calculate average node degree"
            )

        return sum(self.graph.degree()) / self.graph.vcount()

    def calculate_total_length(self) -> int:
        """
        Calculate the total length of all contigs in the graph.

        Returns
        -------
        int
            Total length of all contigs.

        Examples
        --------
        >>> cg.calculate_total_length()
        120000
        """
        contig_lengths = [
            len(self.contig_parser.get_sequence(seq)) for seq in self.contig_names
        ]
        return sum(contig_lengths)

    def calculate_average_contig_length(self) -> float:
        """
        Calculate the average contig length.

        Returns
        -------
        int
            Average contig length.

        Raises
        ------
        ValueError
            If the graph does not have any contig.

        Examples
        --------
        >>> cg.calculate_average_contig_length()
        40000
        """

        contig_lengths = [
            len(self.contig_parser.get_sequence(seq)) for seq in self.contig_names
        ]
        if len(contig_lengths) == 0:
            raise ValueError(
                "Graph does not have any contigs, cannot calculate average contig length"
            )

        return sum(contig_lengths) / len(contig_lengths)

    def calculate_n50_l50(self) -> tuple[int, int]:
        """
        Calculate N50 and L50 for the contigs in the graph.

        Returns
        -------
        tuple of (int, int)
            A tuple containing:
            - N50 : int
                The length N such that 50% of the total length is contained in contigs of length >= N.
            - L50 : int
                The minimum number of contigs whose summed length >= 50% of the total length.

        Examples
        --------
        >>> cg.calculate_n50_l50()
        (15000, 12)
        """

        contig_lengths = [
            len(self.contig_parser.get_sequence(seq)) for seq in self.contig_names
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
        Calculate the GC content of contig sequences.

        Returns
        -------
        float
            GC content as a percentage of total base pairs.

        Raises
        ------
        ValueError
            If total length of the contigs is zero.

        Examples
        --------
        >>> cg.get_gc_content()
        0.42
        """

        contig_sequences = [
            self.contig_parser.get_sequence(seq) for seq in self.contig_names
        ]
        total_length = self.calculate_total_length()

        if total_length == 0:
            raise ValueError(
                "Total length of contigs is zero, cannot calculate GC content"
            )

        gc_count = sum(seq.count("G") + seq.count("C") for seq in contig_sequences)
        return gc_count / total_length
