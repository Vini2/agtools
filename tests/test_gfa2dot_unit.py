#!/usr/bin/env python3

import importlib
import types

from agtools.commands.gfa2dot import _write_abyss_dot, _write_dot, gfa2dot

gfa2dot_module = importlib.import_module("agtools.commands.gfa2dot")


def test_write_abyss_dot_emits_nodes_and_links(tmp_path):
    dummy_graph = types.SimpleNamespace(
        graph=types.SimpleNamespace(vs={"name": ["seg1", "seg2"]}),
        segment_lengths={"seg1": 4, "seg2": 6},
        link_overlap={("seg1+", "seg2-"): 5},
    )

    output_file = _write_abyss_dot(dummy_graph, str(tmp_path))
    content = (tmp_path / "graph.gv").read_text()

    assert output_file == str(tmp_path / "graph.gv")
    assert '"seg1+" [l=4]' in content
    assert '"seg2-" [l=6]' in content
    assert '"seg1+" -> "seg2-" [d=-5]' in content


def test_write_dot_delegates_to_igraph_writer(tmp_path):
    class DummyInnerGraph:
        def __init__(self):
            self.called_path = None

        def write_dot(self, path):
            self.called_path = path
            with open(path, "w") as f:
                f.write("graph {}")

    inner_graph = DummyInnerGraph()
    wrapper = types.SimpleNamespace(graph=inner_graph)

    output_file = _write_dot(wrapper, str(tmp_path))

    assert output_file == str(tmp_path / "graph.dot")
    assert inner_graph.called_path == str(tmp_path / "graph.dot")
    assert (tmp_path / "graph.dot").read_text() == "graph {}"


def test_gfa2dot_selects_writer_based_on_flag(tmp_path, monkeypatch):
    dummy_ug = object()
    calls = []

    monkeypatch.setattr(
        gfa2dot_module.UnitigGraph, "from_gfa", staticmethod(lambda _path: dummy_ug)
    )

    def fake_abyss(graph, output_path):
        calls.append(("abyss", graph, output_path))
        return "abyss-output"

    def fake_dot(graph, output_path):
        calls.append(("dot", graph, output_path))
        return "dot-output"

    monkeypatch.setattr(gfa2dot_module, "_write_abyss_dot", fake_abyss)
    monkeypatch.setattr(gfa2dot_module, "_write_dot", fake_dot)

    assert gfa2dot("graph.gfa", abyss=True, output_path=str(tmp_path)) == "abyss-output"
    assert gfa2dot("graph.gfa", abyss=False, output_path=str(tmp_path)) == "dot-output"
    assert calls == [
        ("abyss", dummy_ug, str(tmp_path)),
        ("dot", dummy_ug, str(tmp_path)),
    ]
