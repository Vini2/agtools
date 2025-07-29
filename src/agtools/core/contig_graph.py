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
