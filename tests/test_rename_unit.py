#!/usr/bin/env python3

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

    output_file = _write_renamed_file(
        str(gfa_file),
        {"seg1": "pref_seg1", "seg2": "pref_seg2"},
        {"path1": "pref_path1"},
        {"walk1": "pref_walk1"},
        str(tmp_path),
    )

    content = (tmp_path / "renamed_graph.gfa").read_text()

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

    output_file = rename(str(gfa_file), "pref", str(tmp_path))

    content = (tmp_path / "renamed_graph.gfa").read_text()
    assert output_file == str(tmp_path / "renamed_graph.gfa")
    assert "S\tpref_seg1\tATGC" in content
    assert "P\tpref_path1\tpref_seg1+\t*" in content

