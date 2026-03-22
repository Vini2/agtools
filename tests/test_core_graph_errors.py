#!/usr/bin/env python3

import pytest
from Bio.Seq import Seq
from igraph import Graph

from agtools.core.contig_graph import ContigGraph
from agtools.core.unitig_graph import UnitigGraph


class DummyParser:
    def __init__(self, sequences=None):
        self.sequences = sequences or {}

    def get_sequence(self, seq_name):
        return Seq(self.sequences.get(seq_name, ""))


def _build_empty_contig_graph():
    return ContigGraph(
        graph=Graph(directed=False),
        vcount=0,
        lcount=0,
        ecount=0,
        file_path="",
        contig_names=[],
        contig_name_to_id={},
        contig_parser=DummyParser(),
        contig_descriptions=None,
        graph_to_contig_map=None,
        self_loops=[],
    )


def test_contig_graph_average_degree_raises_on_empty_graph():
    contig_graph = _build_empty_contig_graph()

    with pytest.raises(ValueError, match="does not have any contigs"):
        contig_graph.calculate_average_node_degree()


def test_contig_graph_average_length_raises_on_empty_graph():
    contig_graph = _build_empty_contig_graph()

    with pytest.raises(ValueError, match="does not have any contigs"):
        contig_graph.calculate_average_contig_length()


def test_contig_graph_gc_content_raises_on_zero_length():
    contig_graph = _build_empty_contig_graph()

    with pytest.raises(ValueError, match="Total length of contigs is zero"):
        contig_graph.get_gc_content()


def test_contig_graph_is_connected_raises_for_unknown_names():
    contig_graph = _build_empty_contig_graph()

    with pytest.raises(KeyError, match="Contig names do not exist in the assembly"):
        contig_graph.is_connected("contig_a", "contig_b")


def test_contig_graph_adjacency_matrix_rejects_invalid_type():
    contig_graph = _build_empty_contig_graph()

    with pytest.raises(ValueError, match="type must be 'matrix' or 'pandas'"):
        contig_graph.get_adjacency_matrix(type="invalid")


def test_unitig_graph_average_degree_raises_on_empty_graph():
    unitig_graph = UnitigGraph()

    with pytest.raises(ValueError, match="does not have any segments"):
        unitig_graph.calculate_average_node_degree()


def test_unitig_graph_average_length_raises_on_empty_graph():
    unitig_graph = UnitigGraph()

    with pytest.raises(ValueError, match="does not have any segments"):
        unitig_graph.calculate_average_segment_length()


def test_unitig_graph_gc_content_raises_on_zero_length():
    unitig_graph = UnitigGraph()

    with pytest.raises(ValueError, match="Total length of segments is zero"):
        unitig_graph.get_gc_content()


def test_unitig_graph_adjacency_matrix_rejects_invalid_type():
    unitig_graph = UnitigGraph()

    with pytest.raises(ValueError, match="type must be 'matrix' or 'pandas'"):
        unitig_graph.get_adjacency_matrix(type="invalid")


def test_unitig_graph_get_path_raises_for_unknown_path():
    unitig_graph = UnitigGraph()

    with pytest.raises(KeyError, match="Unknown path"):
        unitig_graph.get_path("missing_path")
