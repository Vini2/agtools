#!/usr/bin/env python3

import os
import tempfile

from agtools.core.graph import UnitigGraph

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]


def test_from_gfa_basic_segments_and_links():
    gfa_content = "S\tseg1\tATGC\nS\tseg2\tGGTT\nL\tseg1\t+\tseg2\t-\t10M\n"

    with tempfile.NamedTemporaryFile("w+", delete=False) as f:
        f.write(gfa_content)
        f_path = f.name

    ug = UnitigGraph.from_gfa(f_path)
    os.unlink(f_path)

    # Check segments
    assert "seg1" in ug.segment_sequences
    assert "seg2" in ug.segment_sequences
    assert ug.graph.vcount() == len(ug.segment_sequences)

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
