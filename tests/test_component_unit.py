#!/usr/bin/env python3

import importlib
import types

from agtools.commands.component import _write_component_graph, component

component_module = importlib.import_module("agtools.commands.component")


def test_write_component_graph_keeps_only_requested_component(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "H\tVN:Z:1.0\n"
        "S\ts1\tAAAA\n"
        "S\ts2\tTT\n"
        "L\ts1\t+\ts2\t-\t1M\n"
        "C\ts1\t+\ts2\t+\t1M\n"
        "C\ts1\t+\ts1\t+\t1M\n"
        "P\tp1\ts1+,s2-\t*\n"
        "P\tp2\ts1+\t*\n"
        "W\tw1\t*\t*\t*\t>s1<s2\n"
        "W\tw2\t*\t*\t*\t>s1\n"
    )

    target = tmp_path / "component_graph.gfa"
    output_file = _write_component_graph({"s1"}, str(gfa_file), str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "component_graph.gfa")
    assert "S\ts1\tAAAA" in content
    assert "S\ts2\tTT" not in content
    assert "L\ts1\t+\ts2\t-\t1M" not in content
    assert "C\ts1\t+\ts2\t+\t1M" not in content
    assert "C\ts1\t+\ts1\t+\t1M" in content
    assert "P\tp1\ts1+,s2-\t*" not in content
    assert "W\tw1\t*\t*\t*\t>s1<s2" not in content
    assert "P\tp2\ts1+\t*" in content
    assert "W\tw2\t*\t*\t*\t>s1" in content


def test_component_selects_component_containing_requested_segment(tmp_path, monkeypatch):
    class DummyUG:
        def __init__(self):
            self.graph = types.SimpleNamespace(components=lambda: [[0, 2], [1]])
            self.segment_name_to_id = {"s1": 0, "s2": 1, "s3": 2}
            self.segment_names = ["s1", "s2", "s3"]

    monkeypatch.setattr(
        component_module.UnitigGraph, "from_gfa", staticmethod(lambda _: DummyUG())
    )

    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "S\ts1\tAAAA\n"
        "S\ts2\tTT\n"
        "S\ts3\tGG\n"
        "L\ts1\t+\ts3\t+\t1M\n"
        "L\ts2\t+\ts1\t+\t1M\n"
    )

    target = tmp_path / "component_graph.gfa"
    output_file = component(str(gfa_file), segment="s3", output_path=str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "component_graph.gfa")
    assert "S\ts1\tAAAA" in content
    assert "S\ts3\tGG" in content
    assert "S\ts2\tTT" not in content
    assert "L\ts1\t+\ts3\t+\t1M" in content
    assert "L\ts2\t+\ts1\t+\t1M" not in content
