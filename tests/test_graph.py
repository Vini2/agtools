#!/usr/bin/env python3

import os
import tempfile

from agtools.core.graph import UnitigGraph


def test_from_gfa_basic_segments_and_links():
    gfa_content = "S\tseg1\tATGC\nS\tseg2\tGGTT\nL\tseg1\t+\tseg2\t-\t10M\n"

    with tempfile.NamedTemporaryFile("w+", delete=False) as f:
        f.write(gfa_content)
        f_path = f.name

    ag = UnitigGraph.from_gfa(f_path)
    os.unlink(f_path)

    # Check segments
    assert "seg1" in ag.segment_sequences
    assert "seg2" in ag.segment_sequences

    # Check segment lengths
    assert ag.segment_lengths["seg1"] == 4
    assert ag.segment_lengths["seg2"] == 4

    # Check neighbours
    assert ag.get_neighbors("seg1") == ["seg2"]
    assert ag.get_neighbors("seg2") == ["seg1"]

    # Check edge presence
    assert ag.graph.ecount() == 1
    assert ag.graph.vcount() == 2


def test_oriented_links_and_overlap():
    gfa_content = "S\tsegA\tACTG\nS\tsegB\tTGCA\nL\tsegA\t+\tsegB\t-\t5M\n"

    with tempfile.NamedTemporaryFile("w+", delete=False) as f:
        f.write(gfa_content)
        f_path = f.name

    ag = UnitigGraph.from_gfa(f_path)
    os.unlink(f_path)

    # Oriented links should be symmetric
    assert ag.oriented_links["segA"]["segB"] == [("+", "-")]
    assert ag.oriented_links["segB"]["segA"] == [("+", "-")]

    # Overlap stored in both orientations
    assert ag.link_overlap[("segA+", "segB-")] == 5
    assert ag.link_overlap[("segB+", "segA-")] == 5


def test_self_loops_are_recorded():
    gfa_content = "S\tsegX\tGATTACA\nL\tsegX\t+\tsegX\t-\t7M\n"

    with tempfile.NamedTemporaryFile("w+", delete=False) as f:
        f.write(gfa_content)
        f_path = f.name

    ag = UnitigGraph.from_gfa(f_path)
    os.unlink(f_path)

    assert "segX" in ag.self_loops
    assert ag.graph.ecount() == 0  # loop removed by simplify()
    assert ag.graph.vcount() == 1
