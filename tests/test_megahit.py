#!/usr/bin/env python3

import pathlib

import pytest

from agtools.assemblers import megahit

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def contig_graph():
    """Load the contig graph once per test module."""
    graph_file = DATADIR / "5G" / "final.gfa"
    contigs_file = DATADIR / "5G" / "final.contigs.fa"
    return megahit.get_contig_graph(graph_file, contigs_file)


def test_get_contig_graph(contig_graph):

    assert contig_graph.vcount == 11761
    assert contig_graph.ecount == 1120

    assert len(contig_graph.contig_names) == 11761


def test_contig_names_mappings(contig_graph):

    assert contig_graph.contig_names[0] == "NODE_1_length_205_cov_1.0000_ID_1"
    assert (
        contig_graph.graph_to_contig_map["NODE_1_length_205_cov_1.0000_ID_1"]
        == "k141_4704"
    )
    assert (
        contig_graph.contig_descriptions["k141_4704"]
        == "k141_4704 flag=0 multi=1.0000 len=205"
    )


def test_contig_get_neighbors(contig_graph):

    assert "NODE_8191_length_113647_cov_17.9932_ID_16381" in contig_graph.get_neighbors(
        "NODE_9687_length_564_cov_18.0000_ID_19373"
    )
    assert "NODE_8333_length_892_cov_32.0000_ID_16665" in contig_graph.get_neighbors(
        "NODE_11566_length_60808_cov_17.9787_ID_23131"
    )


def test_is_connected(contig_graph):

    assert contig_graph.is_connected(
        "NODE_8333_length_892_cov_32.0000_ID_16665",
        "NODE_3504_length_251_cov_17.7273_ID_7007",
    )
    assert contig_graph.is_connected(
        "NODE_9687_length_564_cov_18.0000_ID_19373",
        "NODE_10530_length_682131_cov_19.0000_ID_21059",
    )


def test_contig_sequences(contig_graph):

    assert contig_graph.contig_parser.get_sequence(
        "NODE_11761_length_301_cov_1.0000_ID_23521"
    ).startswith("GCCGATGCCGCC")
    assert contig_graph.contig_parser.get_sequence(
        "NODE_1_length_205_cov_1.0000_ID_1"
    ).endswith("AAAATGACCCGAA")


def test_contig_index(contig_graph):

    assert contig_graph.contig_parser.index["k141_11282"] == 30816923
    assert contig_graph.contig_parser.index["k141_8112"] == 28859669


def test_adjacency_matrix(contig_graph):

    assert len(contig_graph.get_adjacency_matrix()) == 11761
    assert contig_graph.get_adjacency_matrix()[9686, 8190] == 1


def test_connected_components(contig_graph):
    assert len(contig_graph.get_connected_components()) == 10657


def test_average_node_degree(contig_graph):
    assert contig_graph.calculate_average_node_degree() == 0


def test_total_length(contig_graph):
    assert contig_graph.calculate_total_length() == 30368344


def test_average_contig_length(contig_graph):
    assert contig_graph.calculate_average_contig_length() == 2582


def test_n50_l50(contig_graph):
    assert contig_graph.calculate_n50_l50() == (151640, 52)


def test_gc_content(contig_graph):
    assert contig_graph.get_gc_content() == 0.6378570066250566
