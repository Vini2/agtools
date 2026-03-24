#!/usr/bin/env python3

from agtools.commands.fastg2gfa import _parse_fastg, _write_gfa, fastg2gfa


def test_parse_fastg_extracts_segments_and_oriented_edges(tmp_path):
    fastg_file = tmp_path / "graph.fastg"
    fastg_file.write_text(">1:2,3';\nAAAA\n>2';\nCCCC\n>3;\nGGGG\n")

    segments, edges = _parse_fastg(str(fastg_file))

    assert segments == {"1": "AAAA", "3": "GGGG"}
    assert edges == [("1", "+", "2", "+"), ("1", "+", "3", "-")]


def test_parse_fastg_handles_reverse_from_orientation(tmp_path):
    fastg_file = tmp_path / "graph.fastg"
    fastg_file.write_text(">4':5;\nTT\n>5;\nAA\n")

    segments, edges = _parse_fastg(str(fastg_file))

    assert segments == {"5": "AA"}
    assert edges == [("4", "-", "5", "+")]


def test_write_gfa_outputs_header_segments_and_links(tmp_path):
    segments = {"seg1": "ATGC", "seg2": "GGTT"}
    edges = [("seg1", "+", "seg2", "-",)]

    target = tmp_path / "converted_graph.gfa"
    output_file = _write_gfa(
        segments=segments, edges=edges, output_path=str(target), fixed_overlap=41
    )

    content = target.read_text().splitlines()
    assert output_file == str(target)
    assert content[0] == "H\tVN:Z:1.0"
    assert "S\tseg1\tATGC" in content
    assert "S\tseg2\tGGTT" in content
    assert "L\tseg1\t+\tseg2\t-\t41M" in content


def test_fastg2gfa_end_to_end_generates_expected_gfa(tmp_path):
    fastg_file = tmp_path / "graph.fastg"
    fastg_file.write_text(">A:B';\nACGT\n>B';\nTTTT\n")

    target = tmp_path / "converted_graph.gfa"
    output_file = fastg2gfa(str(fastg_file), k_overlap=55, gfa_path=str(target))

    content = target.read_text().splitlines()

    assert output_file == str(target)
    assert "S\tA\tACGT" in content
    assert "L\tA\t+\tB\t-\t55M" in content
