#!/usr/bin/env python3

import os
import pathlib
import tempfile

import pytest

from agtools.core.unitig_graph import UnitigGraph

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


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

    # Check neighbours
    assert ug.get_neighbours("seg1") == ["seg2"]
    assert ug.get_neighbours("seg2") == ["seg1"]

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
def test_is_connected():
    graph_path = DATADIR / "test_graph.gfa"
    ug = UnitigGraph.from_gfa(graph_path)

    assert not ug.is_connected("seg1", "seg2")
    assert ug.is_connected("seg4", "seg5")
    assert not ug.is_connected("seg10", "segX")


def test_get_sequence_segment():
    graph_path = DATADIR / "test_graph.gfa"
    ug = UnitigGraph.from_gfa(graph_path)

    # Check segment sequences
    assert ug.get_segment_sequence("seg1") == "ATGCGTATGCGTATGCGTAA"
