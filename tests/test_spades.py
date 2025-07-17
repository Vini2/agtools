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

    assert contig_graph.graph.vcount() == 189
    assert contig_graph.graph.ecount() == 394

    assert len(contig_graph.contig_ids) == 189
    assert len(contig_graph.contig_names) == 189

    assert "NODE_1_length_488682_cov_86.190505" in contig_graph.contig_names.values()



def test_get_unitig_graph():

    graph_file = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"

    unitig_graph = spades.get_unitig_graph(graph_file)

    assert unitig_graph.graph.vcount() == 982
    assert unitig_graph.graph.ecount() == 1265