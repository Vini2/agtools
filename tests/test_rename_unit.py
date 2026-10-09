#!/usr/bin/env python3

import pytest

from agtools.commands.rename import (
    _build_element_maps,
    _remap_element,
    _write_renamed_file,
    rename,
)


def test_remap_element_returns_original_if_unmapped():
    assert _remap_element("segX", {"seg1": "pref_seg1"}) == "segX"


def test_build_element_maps_collects_segments_paths_and_walks(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "H\tVN:Z:1.0\n"
        "S\tseg1\tATGC\n"
        "S\tseg2\tGCTA\n"
        "P\tpath1\tseg1+,seg2-\t*\n"
        "W\twalk1\t*\t*\t*\t>seg1<seg2\n"
    )

    segment_map, path_map, walk_map = _build_element_maps(str(gfa_file), "pref")

    assert segment_map == {"seg1": "pref_seg1", "seg2": "pref_seg2"}
    assert path_map == {"path1": "pref_path1"}
    assert walk_map == {"walk1": "pref_walk1"}


def test_write_renamed_file_updates_all_supported_gfa_tags(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "H\tVN:Z:1.0\n"
        "S\tseg1\tATGC\n"
        "S\tseg2\tGCTA\n"
        "L\tseg1\t+\tseg2\t-\t4M\n"
        "J\tseg1\t+\tseg2\t+\t4M\n"
        "C\tseg1\t+\tseg2\t+\t4M\n"
        "P\tpath1\tseg1+,seg2-\t*\n"
        "W\twalk1\t*\t*\t*\t>seg1<seg2\n"
        "X\tcustom\tline\n"
    )

    target = tmp_path / "renamed_graph.gfa"
    output_file = _write_renamed_file(
        str(gfa_file),
        {"seg1": "pref_seg1", "seg2": "pref_seg2"},
        {"path1": "pref_path1"},
        {"walk1": "pref_walk1"},
        str(target),
    )

    content = target.read_text()

    assert output_file == str(tmp_path / "renamed_graph.gfa")
    assert "S\tpref_seg1\tATGC" in content
    assert "L\tpref_seg1\t+\tpref_seg2\t-\t4M" in content
    assert "J\tpref_seg1\t+\tpref_seg2\t+\t4M" in content
    assert "C\tpref_seg1\t+\tpref_seg2\t+\t4M" in content
    assert "P\tpref_path1\tpref_seg1+,pref_seg2-\t*" in content
    assert "W\tpref_walk1\t*\t*\t*\t>pref_seg1<pref_seg2" in content
    assert "X\tcustom\tline" in content


def test_rename_end_to_end(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tATGC\nP\tpath1\tseg1+\t*\n")

    target = tmp_path / "renamed_graph.gfa"
    output_file = rename(str(gfa_file), "pref", str(target))

    content = target.read_text()
    assert output_file == str(tmp_path / "renamed_graph.gfa")
    assert "S\tpref_seg1\tATGC" in content
    assert "P\tpref_path1\tpref_seg1+\t*" in content


def test_build_element_maps_uses_custom_separator(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tATGC\nP\tpath1\tseg1+\t*\nW\twalk1\t*\t*\t*\t>seg1\n")

    segment_map, path_map, walk_map = _build_element_maps(
        str(gfa_file), "pref", separator="."
    )

    assert segment_map == {"seg1": "pref.seg1"}
    assert path_map == {"path1": "pref.path1"}
    assert walk_map == {"walk1": "pref.walk1"}


def test_build_element_maps_defaults_to_underscore(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tATGC\n")

    segment_map, _, _ = _build_element_maps(str(gfa_file), "pref")

    assert segment_map == {"seg1": "pref_seg1"}


def test_build_element_maps_empty_prefix_leaves_ids_unchanged(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tATGC\nP\tpath1\tseg1+\t*\n")

    segment_map, path_map, _ = _build_element_maps(str(gfa_file), "")

    assert segment_map == {"seg1": "seg1"}
    assert path_map == {"path1": "path1"}


def test_build_element_maps_raises_on_collision_with_existing_id(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tATGC\nS\tpref_seg1\tGCTA\n")

    with pytest.raises(ValueError, match="duplicate segment IDs"):
        _build_element_maps(str(gfa_file), "pref")


def test_build_element_maps_collision_depends_on_separator(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tATGC\nS\tpref_seg1\tGCTA\n")

    # A different separator avoids the clash
    segment_map, _, _ = _build_element_maps(str(gfa_file), "pref", separator=".")

    assert segment_map == {"seg1": "pref.seg1", "pref_seg1": "pref.pref_seg1"}


def test_build_element_maps_raises_on_path_collision(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "S\tseg1\tATGC\nP\tpath1\tseg1+\t*\nP\tpref_path1\tseg1+\t*\n"
    )

    with pytest.raises(ValueError, match="duplicate path IDs"):
        _build_element_maps(str(gfa_file), "pref")


def test_build_element_maps_raises_on_walk_collision(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "S\tseg1\tATGC\nW\twalk1\t*\t*\t*\t>seg1\nW\tpref_walk1\t*\t*\t*\t>seg1\n"
    )

    with pytest.raises(ValueError, match="duplicate walk IDs"):
        _build_element_maps(str(gfa_file), "pref")


def test_rename_end_to_end_with_custom_separator(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "S\tseg1\tATGC\n"
        "S\tseg2\tGCTA\n"
        "L\tseg1\t+\tseg2\t-\t4M\n"
        "P\tpath1\tseg1+,seg2-\t*\n"
        "W\twalk1\t*\t*\t*\t>seg1<seg2\n"
    )

    target = tmp_path / "renamed_graph.gfa"
    rename(str(gfa_file), "pref", str(target), separator="-")

    content = target.read_text()
    assert "S\tpref-seg1\tATGC" in content
    assert "L\tpref-seg1\t+\tpref-seg2\t-\t4M" in content
    assert "P\tpref-path1\tpref-seg1+,pref-seg2-\t*" in content
    assert "W\tpref-walk1\t*\t*\t*\t>pref-seg1<pref-seg2" in content


def test_rename_end_to_end_empty_prefix_is_a_no_op(tmp_path):
    original = "S\tseg1\tATGC\nP\tpath1\tseg1+\t*\n"
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(original)

    target = tmp_path / "renamed_graph.gfa"
    rename(str(gfa_file), "", str(target))

    assert target.read_text() == original


def test_rename_raises_on_collision(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tATGC\nS\tpref_seg1\tGCTA\n")

    target = tmp_path / "renamed_graph.gfa"

    with pytest.raises(ValueError, match="duplicate segment IDs"):
        rename(str(gfa_file), "pref", str(target))

    assert not target.exists()
