#!/usr/bin/env python3

import pathlib

from agtools.assemblers import spades

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


def test_get_contig_graph():

    graph_file = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    contig_paths_file = DATADIR / "ESC" / "contigs.paths"

    contig_graph = spades.get_contig_graph(graph_file, contig_paths_file)

    assert contig_graph.vcount == 189
    assert contig_graph.ecount == 394

    assert len(contig_graph.contig_names) == 189

    assert "NODE_1_length_488682_cov_86.190505" in contig_graph.contig_names.values()


def test_contig_get_neighbours():

    graph_file = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    contig_paths_file = DATADIR / "ESC" / "contigs.paths"

    contig_graph = spades.get_contig_graph(graph_file, contig_paths_file)

    assert "NODE_4_length_346431_cov_86.228266" in contig_graph.get_neighbours("NODE_1_length_488682_cov_86.190505")
    assert "NODE_44_length_45842_cov_86.030074" in contig_graph.get_neighbours("NODE_1_length_488682_cov_86.190505")


def test_is_connected():

    graph_file = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    contig_paths_file = DATADIR / "ESC" / "contigs.paths"

    contig_graph = spades.get_contig_graph(graph_file, contig_paths_file)

    assert contig_graph.is_connected("NODE_1_length_488682_cov_86.190505", "NODE_146_length_99_cov_86.818182")
    assert contig_graph.is_connected("NODE_164_length_65_cov_81.100000", "NODE_146_length_99_cov_86.818182")


def test_get_unitig_graph():

    graph_file = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"

    unitig_graph = spades.get_unitig_graph(graph_file)

    assert unitig_graph.vcount == 982
    assert unitig_graph.ecount == 1265
