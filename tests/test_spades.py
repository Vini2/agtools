#!/usr/bin/env python3

import pathlib

import pytest

from agtools.assemblers import spades

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def contig_graph():
    """Load the contig graph once per test module."""
    graph_file = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    contigs_file = DATADIR / "ESC" / "contigs.fasta"
    contig_paths_file = DATADIR / "ESC" / "contigs.paths"
    return spades.get_contig_graph(graph_file, contigs_file, contig_paths_file)


def test_get_contig_graph(contig_graph):

    assert contig_graph.vcount == 189
    assert contig_graph.ecount == 394

    assert len(contig_graph.contig_names) == 189

    assert "NODE_1_length_488682_cov_86.190505" in contig_graph.contig_names.values()


def test_contig_names_mappings(contig_graph):

    assert contig_graph.contig_names[0] == "NODE_1_length_488682_cov_86.190505"
    assert contig_graph.contig_names[1] == "NODE_2_length_472233_cov_17.669606"
    assert contig_graph.contig_names[100] == "NODE_101_length_219_cov_317.097561"


def test_contig_get_neighbors(contig_graph):

    assert "NODE_4_length_346431_cov_86.228266" in contig_graph.get_neighbors(
        "NODE_1_length_488682_cov_86.190505"
    )
    assert "NODE_44_length_45842_cov_86.030074" in contig_graph.get_neighbors(
        "NODE_1_length_488682_cov_86.190505"
    )


def test_is_connected(contig_graph):

    assert contig_graph.is_connected(
        "NODE_1_length_488682_cov_86.190505", "NODE_146_length_99_cov_86.818182"
    )
    assert contig_graph.is_connected(
        "NODE_164_length_65_cov_81.100000", "NODE_146_length_99_cov_86.818182"
    )


def test_contig_sequences(contig_graph):

    assert (
        contig_graph.contig_parser.get_sequence("NODE_174_length_58_cov_650.333333")
        == "GAACTATTATCATTAGCTAAGGTAATAGACAATCAAAGGCTTACCTATTGCTATGCGT"
    )
    assert (
        contig_graph.contig_parser.get_sequence("NODE_189_length_56_cov_33.000000")
        == "TGGCTCTTCAGGATCCAGGGTGTAGTCGGGGTCTGAATCCTCCGGTCTCCAGGAGG"
    )


def test_contig_index(contig_graph):

    assert (
        contig_graph.contig_parser.index["NODE_174_length_58_cov_650.333333"] == 8485847
    )
    assert (
        contig_graph.contig_parser.index["NODE_189_length_56_cov_33.000000"] == 8487228
    )


def test_adjacency_matrix(contig_graph):

    assert len(contig_graph.get_adjacency_matrix()) == 189
    assert contig_graph.get_adjacency_matrix()[3, 0] == 1


def test_connected_components(contig_graph):
    assert len(contig_graph.get_connected_components()) == 3


def test_average_node_degree(contig_graph):
    assert contig_graph.calculate_average_node_degree() == 4


def test_total_length(contig_graph):
    assert contig_graph.calculate_total_length() == 8341464


def test_average_contig_length(contig_graph):
    assert contig_graph.calculate_average_contig_length() == 44134


def test_n50_l50(contig_graph):
    assert contig_graph.calculate_n50_l50() == (220639, 14)


def test_gc_content(contig_graph):
    assert contig_graph.get_gc_content() == 0.4350709899365387


def test_get_unitig_graph():

    graph_file = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"

    unitig_graph = spades.get_unitig_graph(graph_file)

    assert unitig_graph.vcount == 982
    assert unitig_graph.ecount == 1265
