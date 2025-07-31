#!/usr/bin/env python3

import pathlib

import pytest

from agtools.assemblers import flye

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def contig_graph():
    """Load the contig graph once per test module."""
    graph_file = DATADIR / "1Y3B" / "assembly_graph.gfa"
    contig_paths_file = DATADIR / "1Y3B" / "assembly_info.txt"
    contigs_file = DATADIR / "1Y3B" / "assembly.fasta"
    return flye.get_contig_graph(graph_file, contigs_file, contig_paths_file)


def test_get_contig_graph(contig_graph):

    assert contig_graph.vcount == 67
    assert contig_graph.ecount == 2

    assert len(contig_graph.contig_names) == 67

    assert "contig_6" in contig_graph.contig_names.values()


def test_is_connected(contig_graph):

    assert contig_graph.is_connected("contig_64", "contig_9")
    assert contig_graph.is_connected("contig_64", "contig_3") == False


def test_contig_sequences(contig_graph):

    assert contig_graph.contig_parser.get_sequence("contig_2").startswith(
        "GAATTATAATTTGAAA"
    )
    assert contig_graph.contig_parser.get_sequence("contig_57").endswith(
        "ATATCATCTGATG"
    )


def test_contig_index(contig_graph):

    assert contig_graph.contig_parser.index["contig_2"] == 8896084
    assert contig_graph.contig_parser.index["contig_57"] == 20574852


def test_contig_mappings(contig_graph):

    assert contig_graph.contig_names[0] == "contig_1"
    assert contig_graph.contig_names[5] == "contig_6"


def test_adjacency_matrix(contig_graph):

    assert len(contig_graph.get_adjacency_matrix()) == 67
    assert contig_graph.get_adjacency_matrix()[8, 61] == 1


def test_connected_components(contig_graph):
    assert len(contig_graph.get_connected_components()) == 65


def test_average_node_degree(contig_graph):
    assert contig_graph.calculate_average_node_degree() == 0


def test_total_length(contig_graph):
    assert contig_graph.calculate_total_length() == 23184836


def test_average_contig_length(contig_graph):
    assert contig_graph.calculate_average_contig_length() == 346042


def test_n50_l50(contig_graph):
    assert contig_graph.calculate_n50_l50() == (752071, 4)


def test_gc_content(contig_graph):
    assert contig_graph.get_gc_content() == 0.4752605107924852


def test_get_unitig_graph():

    graph_file = DATADIR / "1Y3B" / "assembly_graph.gfa"

    unitig_graph = flye.get_unitig_graph(graph_file)

    assert unitig_graph.vcount == 69
    assert unitig_graph.ecount == 4
