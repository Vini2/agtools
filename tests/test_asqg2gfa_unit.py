#!/usr/bin/env python3

import pytest

from agtools.commands.asqg2gfa import _get_segments_and_links, _write_gfa, asqg2gfa


def test_get_segments_and_links_parses_orientations_and_filters_invalid_overlaps(tmp_path):
    asqg_file = tmp_path / "graph.asqg"
    asqg_file.write_text(
        "VT\tcontig1\tAAAA\n"
        "VT\tcontig2\tCCCC\n"
        "VT\tcontig3\tGGGG\n"
        "ED\tcontig1 contig2 0 10 0 0 10 0 1\n"
        "ED\tcontig2 contig3 5 15 0 3 13 0 0\n"
        "ED\tcontig1 contig3 0 10 0 0 8 0 1\n"
    )

    segments, links = _get_segments_and_links(str(asqg_file))

    assert segments == {"contig1": "AAAA", "contig2": "CCCC", "contig3": "GGGG"}
    assert links == [
        ["contig1", "+", "contig2", "-", 10],
        ["contig2", "+", "contig3", "+", 10],
    ]


def test_write_gfa_outputs_segments_and_links(tmp_path):
    target = tmp_path / "converted_graph.gfa"
    output_file = _write_gfa(
        segments={"c1": "ATGC", "c2": "GGTT"},
        links=[["c1", "+", "c2", "-", 4]],
        output_path=str(target),
    )

    content = target.read_text().splitlines()

    assert output_file == str(target)
    assert "S\tc1\tATGC" in content
    assert "S\tc2\tGGTT" in content
    assert "L\tc1\t+\tc2\t-\t4M" in content


def test_asqg2gfa_end_to_end(tmp_path):
    asqg_file = tmp_path / "graph.asqg"
    asqg_file.write_text(
        "VT\tcontig1\tAAAA\n" "VT\tcontig2\tCCCC\n" "ED\tcontig1 contig2 0 3 0 0 3 0 1\n"
    )

    target = tmp_path / "converted_graph.gfa"
    output_file = asqg2gfa(str(asqg_file), str(target))
    content = target.read_text()

    assert output_file == str(target)
    assert "S\tcontig1\tAAAA" in content
    assert "S\tcontig2\tCCCC" in content
    assert "L\tcontig1\t+\tcontig2\t-\t3M" in content


def test_get_segments_and_links_raises_on_malformed_ed_line(tmp_path):
    asqg_file = tmp_path / "graph.asqg"
    asqg_file.write_text("VT\tcontig1\tAAAA\nED\tcontig1 contig2 0 3\n")

    with pytest.raises(ValueError, match="Malformed ED line"):
        _get_segments_and_links(str(asqg_file))
