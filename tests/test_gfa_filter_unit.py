#!/usr/bin/env python3

from agtools.core.gfa_filter import (
    _ensure_newline,
    parse_path_segment_ids,
    parse_walk_segment_ids,
    write_filtered_gfa,
)


def test_parse_path_segment_ids_strips_orientation_and_supports_multiple_separators():
    assert parse_path_segment_ids("seg1+,seg2-;seg3+") == ["seg1", "seg2", "seg3"]


def test_parse_walk_segment_ids_splits_on_walk_orientation_markers():
    assert parse_walk_segment_ids(">seg1<seg2>seg3") == ["seg1", "seg2", "seg3"]


def test_ensure_newline_only_appends_when_needed():
    assert _ensure_newline("S\tseg1\tAAAA") == "S\tseg1\tAAAA\n"
    assert _ensure_newline("S\tseg1\tAAAA\n") == "S\tseg1\tAAAA\n"


def test_write_filtered_gfa_handles_supported_tags_and_transformations(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "\n"
        "S\tkeep\tAAAA\tLN:i:4\n"
        "S\tdrop\tTT\n"
        "J\tkeep\t+\tdrop\t+\t*\n"
        "C\tkeep\t+\tkeep\t+\t0\t0\t0M\n"
        "P\tpath_drop\tkeep+,drop-\t*\n"
        "P\tpath_keep\tkeep+\t*\n"
        "W\twalk_drop\t*\t*\t*\t>keep<drop\n"
        "W\twalk_keep\t*\t*\t*\t>keep\n"
        "X\tcustom\tfield\n"
    )

    target = tmp_path / "filtered.gfa"

    def keep_segment(seg_id: str) -> bool:
        return seg_id != "drop"

    def transform_segment(parts: list[str]) -> list[str]:
        return [parts[0], parts[1], parts[2].lower(), *parts[3:]]

    write_filtered_gfa(
        str(gfa_file),
        str(target),
        keep_segment=keep_segment,
        transform_segment=transform_segment,
    )

    content = target.read_text()

    assert content.startswith("\n")
    assert "S\tkeep\taaaa\tLN:i:4\n" in content
    assert "S\tdrop\tTT\n" not in content
    assert "J\tkeep\t+\tdrop\t+\t*\n" not in content
    assert "C\tkeep\t+\tkeep\t+\t0\t0\t0M\n" in content
    assert "P\tpath_drop\tkeep+,drop-\t*\n" not in content
    assert "P\tpath_keep\tkeep+\t*\n" in content
    assert "W\twalk_drop\t*\t*\t*\t>keep<drop\n" not in content
    assert "W\twalk_keep\t*\t*\t*\t>keep\n" in content
    assert "X\tcustom\tfield\n" in content
