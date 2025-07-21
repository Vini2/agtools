#!/usr/bin/env python3

import pathlib

from agtools.assemblers import megahit

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


def test_get_contig_graph():

    graph_file = DATADIR / "5G" / "final.gfa"
    contig_file = DATADIR / "5G" / "final.contigs.fa"

    contig_graph = megahit.get_contig_graph(graph_file, contig_file)

    assert contig_graph.vcount == 11761
    assert contig_graph.ecount == 1120

    assert len(contig_graph.contig_sequences) == 11761
    assert len(contig_graph.contig_names) == 11761
