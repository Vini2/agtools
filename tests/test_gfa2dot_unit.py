#!/usr/bin/env python3

import importlib
import types

import pytest

from agtools.commands.gfa2dot import _write_abyss_dot, _write_dot, gfa2dot

gfa2dot_module = importlib.import_module("agtools.commands.gfa2dot")


def test_write_abyss_dot_emits_nodes_and_links(tmp_path):
    dummy_graph = types.SimpleNamespace(
        graph=types.SimpleNamespace(vs={"name": ["seg1", "seg2"]}),
        segment_names=["seg1", "seg2"],
        segment_lengths={"seg1": 4, "seg2": 6},
        link_overlap={(0, "+", 1, "-"): 5, (1, "+", 0, "-"): 5},
    )

    target = tmp_path / "graph.gv"
    output_file = _write_abyss_dot(dummy_graph, str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "graph.gv")
    assert '"seg1+" [l=4]' in content
    assert '"seg2-" [l=6]' in content
    assert '"seg1+" -> "seg2-" [d=-5]' in content
    assert '"seg2+" -> "seg1-" [d=-5]' in content


def test_write_dot_delegates_to_igraph_writer(tmp_path):
    class DummyInnerGraph:
        def __init__(self):
            self.called_handle = None

        def write_dot(self, handle):
            self.called_handle = handle
            handle.write("graph {}")

    inner_graph = DummyInnerGraph()
    wrapper = types.SimpleNamespace(graph=inner_graph)

    target = tmp_path / "graph.dot"
    output_file = _write_dot(wrapper, str(target))

    assert output_file == str(target)
    assert inner_graph.called_handle is not None
    assert target.read_text() == "graph {}"


def test_write_dot_supports_stdout(capsys):
    class DummyInnerGraph:
        def write_dot(self, handle):
            handle.write("graph {}")

    wrapper = types.SimpleNamespace(graph=DummyInnerGraph())

    output_file = _write_dot(wrapper, "-")
    stdout = capsys.readouterr().out

    assert output_file == "-"
    assert stdout == "graph {}"


def test_gfa2dot_selects_writer_based_on_flag(tmp_path, monkeypatch):
    dummy_ug = object()
    calls = []

    monkeypatch.setattr(
        gfa2dot_module.UnitigGraph, "from_gfa", staticmethod(lambda _path: dummy_ug)
    )
    monkeypatch.setattr(gfa2dot_module, "validate_gfa_input", lambda *_args: None)

    def fake_abyss(graph, output_path):
        calls.append(("abyss", graph, output_path))
        return "abyss-output"

    def fake_dot(graph, output_path):
        calls.append(("dot", graph, output_path))
        return "dot-output"

    monkeypatch.setattr(gfa2dot_module, "_write_abyss_dot", fake_abyss)
    monkeypatch.setattr(gfa2dot_module, "_write_dot", fake_dot)

    abyss_target = tmp_path / "graph.gv"
    dot_target = tmp_path / "graph.dot"
    assert (
        gfa2dot("graph.gfa", abyss=True, output_path=str(abyss_target))
        == "abyss-output"
    )
    assert (
        gfa2dot("graph.gfa", abyss=False, output_path=str(dot_target)) == "dot-output"
    )
    assert calls == [
        ("abyss", dummy_ug, str(abyss_target)),
        ("dot", dummy_ug, str(dot_target)),
    ]


def test_gfa2dot_rejects_fastg_input(tmp_path):
    fastg_file = tmp_path / "graph.fastg"
    fastg_file.write_text(">A:B';\nACGT\n>B';\nTTTT\n")

    with pytest.raises(ValueError, match="looks like a FASTG file"):
        gfa2dot(str(fastg_file), abyss=False, output_path=str(tmp_path / "graph.dot"))
