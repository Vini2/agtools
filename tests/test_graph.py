#!/usr/bin/env python3

import os
import pathlib
import tempfile

import pytest

from agtools.core.unitig_graph import UnitigGraph

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def unitig_graph():
    """Load the contig graph once per test module."""
    graph_file = DATADIR / "test_graph.gfa"
    return UnitigGraph.from_gfa(graph_file)


@pytest.fixture(scope="module")
def spades_unitig_graph():
    """Load the contig graph once per test module."""
    graph_file = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    return UnitigGraph.from_gfa(graph_file)


def test_from_gfa_basic_segments_and_links():
    gfa_content = "S\tseg1\tATGC\nS\tseg2\tGGTT\nL\tseg1\t+\tseg2\t-\t10M\n"

    with tempfile.NamedTemporaryFile("w+", delete=False) as f:
        f.write(gfa_content)
        f_path = f.name

    ug = UnitigGraph.from_gfa(f_path)
    os.unlink(f_path)

    # Check segments
    assert "seg1" in ug.segment_lengths
    assert "seg2" in ug.segment_lengths
    assert ug.graph.vcount() == len(ug.segment_lengths)

    # Check segment lengths
    assert ug.segment_lengths["seg1"] == 4
    assert ug.segment_lengths["seg2"] == 4

    # Check neighbors
    assert ug.get_neighbors("seg1") == ["seg2"]
    assert ug.get_neighbors("seg2") == ["seg1"]

    # Check edge presence
    assert ug.graph.ecount() == 1
    assert ug.graph.vcount() == 2


def test_oriented_links_and_overlap():
    gfa_content = "S\tsegA\tACTG\nS\tsegB\tTGCA\nL\tsegA\t+\tsegB\t-\t5M\n"

    with tempfile.NamedTemporaryFile("w+", delete=False) as f:
        f.write(gfa_content)
        f_path = f.name

    ug = UnitigGraph.from_gfa(f_path)
    os.unlink(f_path)

    # Oriented links should be symmetric
    assert ug.oriented_links["segA"]["segB"] == [("+", "-")]
    assert ug.oriented_links["segB"]["segA"] == [("+", "-")]

    # Overlap stored in both orientations
    assert ug.link_overlap[("segA+", "segB-")] == 5
    assert ug.link_overlap[("segB+", "segA-")] == 5


def test_self_loops_are_recorded():
    gfa_content = "S\tsegX\tGATTACA\nL\tsegX\t+\tsegX\t-\t7M\n"

    with tempfile.NamedTemporaryFile("w+", delete=False) as f:
        f.write(gfa_content)
        f_path = f.name

    ug = UnitigGraph.from_gfa(f_path)
    os.unlink(f_path)

    assert "segX" in ug.self_loops
    assert ug.graph.ecount() == 0  # loop removed by simplify()
    assert ug.graph.vcount() == 1


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_is_connected(unitig_graph):
    assert not unitig_graph.is_connected("seg1", "seg2")
    assert unitig_graph.is_connected("seg4", "seg5")
    assert not unitig_graph.is_connected("seg10", "segX")


def test_get_sequence_segment(unitig_graph):
    assert unitig_graph.get_segment_sequence("seg1") == "ATGCGTATGCGTATGCGTAA"


def test_graph_stats(spades_unitig_graph):
    assert spades_unitig_graph.vcount == 982
    assert spades_unitig_graph.ecount == 1265
    assert len(spades_unitig_graph.self_loops) == 1


def test_connected_components(spades_unitig_graph):
    assert len(spades_unitig_graph.get_connected_components()) == 4


def test_average_node_degree(spades_unitig_graph):
    assert spades_unitig_graph.calculate_average_node_degree() == 2


def test_total_length(spades_unitig_graph):
    assert spades_unitig_graph.calculate_total_length() == 8337494


def test_average_segment_length(spades_unitig_graph):
    assert spades_unitig_graph.calculate_average_segment_length() == 8490


def test_n50_l50(spades_unitig_graph):
    assert spades_unitig_graph.calculate_n50_l50() == (60706, 36)


def test_gc_content(spades_unitig_graph):
    assert spades_unitig_graph.get_gc_content() == 0.4351642112126258
