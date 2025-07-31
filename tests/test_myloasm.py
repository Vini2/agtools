#!/usr/bin/env python3

import pathlib

import pytest

from agtools.assemblers import myloasm

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def contig_graph():
    """Load the contig graph once per test module."""
    graph_file = DATADIR / "myloasm" / "final_contig_graph.gfa"
    contigs_file = DATADIR / "myloasm" / "assembly_primary.fa"
    return myloasm.get_contig_graph(graph_file, contigs_file)


def test_get_contig_graph(contig_graph):

    assert contig_graph.vcount == 3
    assert contig_graph.ecount == 1

    assert len(contig_graph.contig_names) == 3


def test_contig_names_mappings(contig_graph):

    assert contig_graph.contig_names[0] == "u673671ctg"
    assert contig_graph.contig_names[1] == "u913838ctg"
    assert contig_graph.contig_names[2] == "u579439ctg"


def test_contig_get_neighbors(contig_graph):

    assert "u579439ctg" in contig_graph.get_neighbors("u913838ctg")


def test_is_connected(contig_graph):

    assert contig_graph.is_connected("u579439ctg", "u913838ctg")


def test_contig_sequences(contig_graph):

    assert contig_graph.contig_parser.get_sequence("u913838ctg").startswith(
        "ATCCTTGCGCATTTTC"
    )
    assert contig_graph.contig_parser.get_sequence("u579439ctg").endswith("GGCGTCGGT")


def test_contig_index(contig_graph):

    assert contig_graph.contig_parser.index["u579439ctg"] == 0
    assert contig_graph.contig_parser.index["u913838ctg"] == 39374


def test_adjacency_matrix(contig_graph):

    assert len(contig_graph.get_adjacency_matrix()) == 3
    assert contig_graph.get_adjacency_matrix()[1, 2] == 1


def test_connected_components(contig_graph):
    assert len(contig_graph.get_connected_components()) == 2


def test_average_node_degree(contig_graph):
    assert contig_graph.calculate_average_node_degree() == 0


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_total_length(contig_graph):
    assert contig_graph.calculate_total_length() == 58400


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_average_contig_length(contig_graph):
    assert contig_graph.calculate_average_contig_length() == 19466


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_n50_l50(contig_graph):
    assert contig_graph.calculate_n50_l50() == (39317, 1)


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_gc_content(contig_graph):
    assert contig_graph.get_gc_content() == 0.551986301369863
