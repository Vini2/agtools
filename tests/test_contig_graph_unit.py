#!/usr/bin/env python3

import pytest
from Bio.Seq import Seq
from igraph import Graph

from agtools.core.contig_graph import ContigGraph


class DummyParser:
    def __init__(self, sequences):
        self.sequences = sequences

    def get_sequence(self, seq_name):
        return Seq(self.sequences[seq_name])


def _build_contig_graph():
    graph = Graph(directed=False)
    graph.add_vertices(3)
    graph.add_edges([(0, 1)])

    contig_names = ["c1", "c2", "c3"]
    contig_name_to_id = {name: i for i, name in enumerate(contig_names)}
    parser = DummyParser({"c1": "GGAA", "c2": "CC", "c3": "AT"})

    return ContigGraph(
        graph=graph,
        vcount=graph.vcount(),
        lcount=1,
        ecount=graph.ecount(),
        file_path="graph.gfa",
        contig_names=contig_names,
        contig_name_to_id=contig_name_to_id,
        contig_parser=parser,
        contig_descriptions=None,
        graph_to_contig_map=None,
        self_loops=[],
    )


def test_contig_graph_query_apis():
    cg = _build_contig_graph()

    assert str(cg.get_contig_sequence("c1")) == "GGAA"
    assert cg.get_neighbors("c1") == ["c2"]
    assert cg.get_neighbors("c3") == []
    assert cg.is_connected("c1", "c2")
    assert not cg.is_connected("c1", "c3")


def test_contig_graph_adjacency_matrix_formats():
    cg = _build_contig_graph()

    adj = cg.get_adjacency_matrix()
    assert adj[0, 1] == 1
    assert adj[0, 2] == 0

    adj_df = cg.get_adjacency_matrix(type="pandas")
    assert list(adj_df.index) == ["c1", "c2", "c3"]
    assert list(adj_df.columns) == ["c1", "c2", "c3"]
    assert adj_df.loc["c1", "c2"] == 1
    assert adj_df.loc["c1", "c3"] == 0


def test_contig_graph_connected_components_and_metrics():
    cg = _build_contig_graph()

    assert len(cg.get_connected_components()) == 2
    assert cg.calculate_average_node_degree() == pytest.approx((1 + 1 + 0) / 3)
    assert cg.calculate_total_length() == 8
    assert cg.calculate_average_contig_length() == pytest.approx(8 / 3)
    assert cg.calculate_n50_l50() == (4, 1)
    assert cg.get_gc_content() == pytest.approx(4 / 8)

