#!/usr/bin/env python3

import importlib
import types

from agtools.commands.stats import _write_stats_file, stats

stats_module = importlib.import_module("agtools.commands.stats")


def test_write_stats_file_contains_expected_metrics(tmp_path):
    stats_dict = {
        "nsegments": 3,
        "nlinks": 2,
        "nloops": 1,
        "ncomponents": 2,
        "average_node_degree": 1.3333333,
        "total_length": 1234,
        "average_segment_length": 411.33,
        "n50": 500,
        "l50": 2,
        "gc_content": 0.4321,
    }

    output_file = _write_stats_file("graph.gfa", stats_dict, str(tmp_path))
    content = (tmp_path / "graph_stats.txt").read_text()

    assert output_file == str(tmp_path / "graph_stats.txt")
    assert "Basic graph statistics for graph.gfa:" in content
    assert "Number of segments: 3" in content
    assert "Number of links: 2" in content
    assert "Number of self-loops: 1" in content
    assert "Sequence-based statistics for graph.gfa:" in content
    assert "N50: 500 bp" in content
    assert "GC content: 43.21%" in content


def test_stats_builds_metrics_from_unitig_graph(tmp_path, monkeypatch):
    class DummyUG:
        def __init__(self):
            self.graph = types.SimpleNamespace(vcount=lambda: 4, ecount=lambda: 3)
            self.self_loops = [0]

        def get_connected_components(self):
            return [[0, 1], [2, 3]]

        def calculate_average_node_degree(self):
            return 1.5

        def calculate_total_length(self):
            return 200

        def calculate_average_segment_length(self):
            return 50

        def calculate_n50_l50(self):
            return (80, 2)

        def get_gc_content(self):
            return 0.5

    monkeypatch.setattr(
        stats_module.UnitigGraph, "from_gfa", staticmethod(lambda _: DummyUG())
    )

    output_file = stats("input.gfa", str(tmp_path))
    content = (tmp_path / "graph_stats.txt").read_text()

    assert output_file == str(tmp_path / "graph_stats.txt")
    assert "Number of segments: 4" in content
    assert "Number of links: 3" in content
    assert "Average node degree: 1.5" in content
    assert "GC content: 50.00%" in content
