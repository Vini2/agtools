#!/usr/bin/env python3

import pytest

from agtools.core.unitig_graph import UnitigGraph


def _write_small_graph(path):
    path.write_text(
        "H\tVN:Z:1.0\n"
        "S\ta\tGGGGAAAA\n"
        "S\tb\tCCAA\n"
        "S\tc\tAT\n"
        "L\ta\t+\tb\t-\t2M\n"
        "L\ta\t+\tb\t-\t2M\n"
        "L\tb\t+\tb\t-\t1M\n"
        "P\tp1\ta+,b-\t*\n"
    )


def test_from_gfa_tracks_paths_and_simplifies_duplicate_links(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    _write_small_graph(gfa_file)

    ug = UnitigGraph.from_gfa(str(gfa_file))

    assert ug.vcount == 3
    assert ug.lcount == 3  # raw link lines in file
    assert ug.ecount == 1  # duplicate edge simplified + self-loop removed
    assert ug.pcount == 1
    assert ug.self_loops == [ug.segment_name_to_id["b"]]


@pytest.mark.filterwarnings("ignore::RuntimeWarning")
def test_unitig_graph_query_and_matrix_apis(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    _write_small_graph(gfa_file)

    ug = UnitigGraph.from_gfa(str(gfa_file))

    assert ug.get_neighbors("a") == ["b"]
    assert ug.get_neighbors("c") == []
    assert ug.is_connected("a", "b")
    assert not ug.is_connected("a", "c")
    assert ug.get_path("p1") == ("a+,b-", "*")

    adj_matrix = ug.get_adjacency_matrix()
    assert adj_matrix[0, 1] == 1
    assert adj_matrix[0, 2] == 0

    adj_df = ug.get_adjacency_matrix(type="pandas")
    assert list(adj_df.index) == ["a", "b", "c"]
    assert list(adj_df.columns) == ["a", "b", "c"]
    assert adj_df.loc["a", "b"] == 1
    assert adj_df.loc["a", "c"] == 0


def test_unitig_graph_sequence_length_and_gc_metrics(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    _write_small_graph(gfa_file)

    ug = UnitigGraph.from_gfa(str(gfa_file))

    assert ug.calculate_total_length() == 14
    assert ug.calculate_average_segment_length() == pytest.approx(14 / 3)
    assert ug.calculate_n50_l50() == (8, 1)
    assert ug.get_gc_content() == pytest.approx(6 / 14)


def test_get_segment_sequence_detects_length_mismatch_after_file_change(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tAAAA\n")

    ug = UnitigGraph.from_gfa(str(gfa_file))

    # Corrupt the backing file after indexing to force a sequence-length mismatch.
    gfa_file.write_text("S\tseg1\tAA\n")

    with pytest.raises(ValueError, match="Wrong sequence retrieved"):
        ug.get_segment_sequence("seg1")
