#!/usr/bin/env python3

import pathlib

from agtools.assemblers import flye

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


def test_get_contig_graph():

    graph_file = DATADIR / "1Y3B" / "assembly_graph.gfa"
    contig_paths_file = DATADIR / "1Y3B" / "assembly_info.txt"
    contigs_file = DATADIR / "1Y3B" / "assembly.fasta"

    contig_graph = flye.get_contig_graph(graph_file, contigs_file, contig_paths_file)

    assert contig_graph.vcount == 67
    assert contig_graph.ecount == 2

    assert len(contig_graph.contig_names) == 67

    assert "contig_6" in contig_graph.contig_names.values()


def test_get_unitig_graph():

    graph_file = DATADIR / "1Y3B" / "assembly_graph.gfa"

    unitig_graph = flye.get_unitig_graph(graph_file)

    assert unitig_graph.vcount == 69
    assert unitig_graph.ecount == 4