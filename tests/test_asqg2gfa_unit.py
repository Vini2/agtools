#!/usr/bin/env python3

import pytest

from agtools.commands.asqg2gfa import _get_segments_and_links, _write_gfa, asqg2gfa
from agtools.commands.gfa2asqg import gfa2asqg


def test_get_segments_and_links_parses_spec_compliant_orientations(tmp_path):
    asqg_file = tmp_path / "graph.asqg"
    asqg_file.write_text(
        "VT\tcontig1\tAAAAAAAA\n"
        "VT\tcontig2\tCCCCCC\n"
        "ED\tcontig1 contig2 5 7 8 0 2 6 0 0\n"
        "ED\tcontig1 contig2 0 1 8 0 1 6 0 0\n"
        "ED\tcontig1 contig2 4 7 8 2 5 6 1 0\n"
        "ED\tcontig1 contig2 0 0 8 5 5 6 1 0\n"
    )

    segments, links = _get_segments_and_links(str(asqg_file))

    assert segments == {"contig1": "AAAAAAAA", "contig2": "CCCCCC"}
    assert links == [
        ["contig1", "+", "contig2", "+", 3],
        ["contig1", "-", "contig2", "+", 2],
        ["contig1", "+", "contig2", "-", 4],
        ["contig1", "-", "contig2", "-", 1],
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
        "VT\tcontig1\tAAAA\n"
        "VT\tcontig2\tCCCC\n"
        "ED\tcontig1 contig2 1 3 4 1 3 4 1 0\n"
    )

    target = tmp_path / "converted_graph.gfa"
    output_file = asqg2gfa(str(asqg_file), str(target))
    content = target.read_text()

    assert output_file == str(target)
    assert "S\tcontig1\tAAAA" in content
    assert "S\tcontig2\tCCCC" in content
    assert "L\tcontig1\t+\tcontig2\t-\t3M" in content


def test_asqg2gfa_round_trips_spec_compliant_gfa2asqg_output(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tA\tAAAA\nS\tB\tCCCC\nL\tA\t+\tB\t-\t3M\n")

    asqg_target = tmp_path / "converted_graph.asqg"
    roundtrip_target = tmp_path / "roundtrip_graph.gfa"

    gfa2asqg(str(gfa_file), str(asqg_target))
    asqg2gfa(str(asqg_target), str(roundtrip_target))

    assert "L\tA\t+\tB\t-\t3M" in roundtrip_target.read_text()


def test_get_segments_and_links_raises_on_malformed_ed_line(tmp_path):
    asqg_file = tmp_path / "graph.asqg"
    asqg_file.write_text("VT\tcontig1\tAAAA\nED\tcontig1 contig2 0 3\n")

    with pytest.raises(ValueError, match="Malformed ED line"):
        _get_segments_and_links(str(asqg_file))


@pytest.mark.parametrize(
    "ed_line",
    [
        "ED\n",
        "ED\tcontig1 contig2 x 3 4 0 3 4 1 0\n",
        "ED\tcontig1 contig2 1 3 4 1 3 4 2 0\n",
        "ED\tcontig1 contig2 1 3 4 0 2 4 1 0\n",
        "ED\tcontig1 contig2 1 3 5 1 3 4 1 0\n",
    ],
)
def test_get_segments_and_links_rejects_other_malformed_ed_variants(tmp_path, ed_line):
    asqg_file = tmp_path / "graph.asqg"
    asqg_file.write_text("VT\tcontig1\tAAAA\nVT\tcontig2\tCCCC\n" + ed_line)

    with pytest.raises(ValueError, match="Malformed ED line"):
        _get_segments_and_links(str(asqg_file))


def test_asqg2gfa_rejects_gfa_input(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tATGC\n")

    with pytest.raises(ValueError, match="looks like a GFA file"):
        asqg2gfa(str(gfa_file), str(tmp_path / "converted_graph.gfa"))
