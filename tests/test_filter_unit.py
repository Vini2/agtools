#!/usr/bin/env python3

import importlib

from agtools.commands.filter import _write_filtered_graph, filter

filter_module = importlib.import_module("agtools.commands.filter")


def test_write_filtered_graph_removes_segments_and_dependent_records(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "H\tVN:Z:1.0\n"
        "S\ts1\tAAAA\n"
        "S\ts2\tTT\n"
        "L\ts1\t+\ts2\t-\t1M\n"
        "J\ts1\t+\ts2\t+\t1M\n"
        "C\ts1\t+\ts2\t+\t1M\n"
        "C\ts1\t+\ts1\t+\t1M\n"
        "P\tp1\ts1+,s2-\t*\n"
        "P\tp2\ts1+\t*\n"
        "W\tw1\t*\t*\t*\t>s1<s2\n"
        "W\tw2\t*\t*\t*\t>s1\n"
    )

    target = tmp_path / "filtered_graph.gfa"
    output_file = _write_filtered_graph({"s2"}, str(gfa_file), str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "filtered_graph.gfa")
    assert "S\ts1\tAAAA" in content
    assert "S\ts2\tTT" not in content
    assert "L\ts1\t+\ts2\t-\t1M" not in content
    assert "J\ts1\t+\ts2\t+\t1M" not in content
    assert "C\ts1\t+\ts2\t+\t1M" not in content
    assert "C\ts1\t+\ts1\t+\t1M" in content
    assert "P\tp1\ts1+,s2-\t*" not in content
    assert "W\tw1\t*\t*\t*\t>s1<s2" not in content
    assert "P\tp2\ts1+\t*" in content
    assert "W\tw2\t*\t*\t*\t>s1" in content


def test_filter_removes_segments_shorter_than_threshold(tmp_path, monkeypatch):
    class DummyUG:
        segment_lengths = {"s1": 10, "s2": 2}

    monkeypatch.setattr(
        filter_module.UnitigGraph, "from_gfa", staticmethod(lambda _: DummyUG())
    )

    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\ts1\tAAAAAAAAAA\nS\ts2\tTT\nP\tp1\ts1+,s2-\t*\n")

    target = tmp_path / "filtered_graph.gfa"
    output_file = filter(str(gfa_file), min_length=5, output_path=str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "filtered_graph.gfa")
    assert "S\ts1\tAAAAAAAAAA" in content
    assert "S\ts2\tTT" not in content
    assert "P\tp1\ts1+,s2-\t*" not in content
