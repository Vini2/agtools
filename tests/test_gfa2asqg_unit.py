#!/usr/bin/env python3

import pytest

from agtools.commands.asqg2gfa import asqg2gfa
from agtools.commands.gfa2asqg import (
    _get_segment_length,
    _get_segments_and_edges,
    _write_asqg,
    gfa2asqg,
)


def test_get_segments_and_edges_converts_link_orientations_to_asqg_coordinates(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "S\tA\tAAAAAAAA\n"
        "S\tB\tCCCCCC\n"
        "L\tA\t+\tB\t+\t3M\n"
        "L\tA\t-\tB\t+\t2M\n"
        "L\tA\t+\tB\t-\t4M\n"
        "L\tA\t-\tB\t-\t1M\n"
    )

    segments, edges = _get_segments_and_edges(str(gfa_file))

    assert segments == {"A": "AAAAAAAA", "B": "CCCCCC"}
    assert edges == [
        ["A", "B", 5, 7, 8, 0, 2, 6, 0, 0],
        ["A", "B", 0, 1, 8, 0, 1, 6, 0, 0],
        ["A", "B", 4, 7, 8, 2, 5, 6, 1, 0],
        ["A", "B", 0, 0, 8, 5, 5, 6, 1, 0],
    ]


def test_write_asqg_outputs_header_segments_and_edges(tmp_path):
    target = tmp_path / "converted_graph.asqg"
    output_file = _write_asqg(
        segments={"A": "ATGC", "B": "GGTT"},
        edges=[["A", "B", 0, 3, 4, 0, 3, 4, 1, 0]],
        output_path=str(target),
    )

    content = target.read_text().splitlines()

    assert output_file == str(target)
    assert content[0] == "HT\tVN:i:1\tER:f:0\tOL:i:4\tCN:i:0\tTE:i:0"
    assert "VT\tA\tATGC" in content
    assert "VT\tB\tGGTT" in content
    assert "ED\tA B 0 3 4 0 3 4 1 0" in content


def test_gfa2asqg_end_to_end_writes_spec_compliant_ed_coordinates(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tA\tAAAA\nS\tB\tCCCC\nL\tA\t+\tB\t-\t3M\n")

    asqg_target = tmp_path / "converted_graph.asqg"
    output_file = gfa2asqg(str(gfa_file), str(asqg_target))
    roundtrip_target = tmp_path / "roundtrip_graph.gfa"
    roundtrip_file = asqg2gfa(str(asqg_target), str(roundtrip_target))

    assert output_file == str(asqg_target)
    assert roundtrip_file == str(roundtrip_target)
    assert "VT\tA\tAAAA" in asqg_target.read_text()
    assert "HT\tVN:i:1\tER:f:0\tOL:i:3\tCN:i:0\tTE:i:0" in asqg_target.read_text()
    assert "ED\tA B 1 3 4 1 3 4 1 0" in asqg_target.read_text()
    assert "L\tA\t+\tB\t-\t3M" in roundtrip_target.read_text()


def test_get_segments_and_edges_rejects_unsupported_overlap_fields(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tA\tAAAA\nS\tB\tCCCC\nL\tA\t+\tB\t+\t*\n")

    with pytest.raises(ValueError, match="Unsupported L overlap field"):
        _get_segments_and_edges(str(gfa_file))


def test_gfa2asqg_rejects_asqg_input(tmp_path):
    asqg_file = tmp_path / "graph.asqg"
    asqg_file.write_text("VT\tseg1\tAAAA\n")

    with pytest.raises(ValueError, match="looks like an ASQG file"):
        gfa2asqg(str(asqg_file), str(tmp_path / "converted_graph.asqg"))


def test_get_segment_length_reads_ln_tag_when_sequence_is_missing():
    parts = ["S", "seg1", "*", "LN:i:12"]

    assert _get_segment_length(parts, "S\tseg1\t*\tLN:i:12\n") == 12


def test_get_segment_length_rejects_invalid_ln_tag():
    parts = ["S", "seg1", "*", "LN:i:not_an_int"]

    with pytest.raises(ValueError, match="Malformed S line"):
        _get_segment_length(parts, "S\tseg1\t*\tLN:i:not_an_int\n")


def test_get_segment_length_rejects_missing_ln_tag_for_sequence_placeholder():
    parts = ["S", "seg1", "*"]

    with pytest.raises(ValueError, match="without sequence data or an LN tag"):
        _get_segment_length(parts, "S\tseg1\t*\n")


def test_get_segments_and_edges_rejects_malformed_s_line(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\n")

    with pytest.raises(ValueError, match="Malformed S line"):
        _get_segments_and_edges(str(gfa_file))


def test_get_segments_and_edges_rejects_malformed_l_line_with_missing_fields(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tA\tAAAA\nL\tA\t+\tB\t+\n")

    with pytest.raises(ValueError, match="Malformed L line"):
        _get_segments_and_edges(str(gfa_file))


def test_get_segments_and_edges_rejects_l_line_with_invalid_orientation(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tA\tAAAA\nS\tB\tCCCC\nL\tA\t?\tB\t+\t1M\n")

    with pytest.raises(ValueError, match="Malformed L line"):
        _get_segments_and_edges(str(gfa_file))


def test_get_segments_and_edges_rejects_links_to_undefined_segments(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tA\tAAAA\nL\tA\t+\tB\t+\t1M\n")

    with pytest.raises(ValueError, match="Link references undefined segment"):
        _get_segments_and_edges(str(gfa_file))


def test_get_segments_and_edges_rejects_overlaps_longer_than_segment_length(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tA\tAAAA\nS\tB\tCCCC\nL\tA\t+\tB\t+\t5M\n")

    with pytest.raises(ValueError, match="Link overlap 5 exceeds segment length"):
        _get_segments_and_edges(str(gfa_file))
