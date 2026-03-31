#!/usr/bin/env python3

import pytest

from agtools.commands.gfa2asqg import _get_segments_and_edges, _write_asqg, gfa2asqg


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

    assert output_file == str(asqg_target)
    assert "VT\tA\tAAAA" in asqg_target.read_text()
    assert "HT\tVN:i:1\tER:f:0\tOL:i:3\tCN:i:0\tTE:i:0" in asqg_target.read_text()
    assert "ED\tA B 1 3 4 1 3 4 1 0" in asqg_target.read_text()


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
