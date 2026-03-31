#!/usr/bin/env python3

import importlib

import pytest

from agtools.commands.fastg2gfa import _parse_fastg
from agtools.commands.gfa2fastg import (
    _get_graph_sequences,
    _write_fastg,
    gfa2fastg,
    reverse_complement,
    reverse_orientation,
)

gfa2fastg_module = importlib.import_module("agtools.commands.gfa2fastg")


def test_reverse_orientation_flips_link_direction():
    assert reverse_orientation("+") == "-"
    assert reverse_orientation("-") == "+"


def test_reverse_complement_returns_expected_sequence():
    assert reverse_complement("AGTC") == "GACT"


def test_get_graph_sequences_extracts_sequences_and_reverse_links(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tA\tAGTC\nS\tB\tTTAA\nL\tA\t+\tB\t-\t55M\n")

    graph_nodes, sequences = _get_graph_sequences(str(gfa_file))

    assert sequences == {"A": "AGTC", "B": "TTAA"}
    assert graph_nodes["A+"] == {"B-"}
    assert graph_nodes["B+"] == {"A-"}


def test_write_fastg_outputs_fastg_headers_and_reverse_complements(tmp_path):
    graph_nodes = {"A+": {"B-"}, "B+": {"A-"}}
    sequences = {"A": "AGTC", "B": "AACC"}

    target = tmp_path / "converted_graph.fastg"
    output_file = _write_fastg(graph_nodes, sequences, str(target))

    content = target.read_text().splitlines()

    assert output_file == str(target)
    assert content == [
        ">A:B';",
        "AGTC",
        ">A';",
        "GACT",
        ">B:A';",
        "AACC",
        ">B';",
        "GGTT",
    ]


def test_gfa2fastg_end_to_end_outputs_parser_compatible_fastg(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tA\tACGT\nS\tB\tTTTT\nL\tA\t+\tB\t-\t55M\n")

    target = tmp_path / "converted_graph.fastg"
    output_file = gfa2fastg(str(gfa_file), str(target))

    segments, edges = _parse_fastg(str(target))

    assert output_file == str(target)
    assert segments == {"A": "ACGT", "B": "TTTT"}
    assert ("A", "+", "B", "-") in edges
    assert ("B", "+", "A", "-") in edges


def test_gfa2fastg_rejects_fastg_input_with_clear_error(tmp_path):
    fastg_file = tmp_path / "graph.fastg"
    fastg_file.write_text(">A:B';\nACGT\n>B';\nTTTT\n")

    with pytest.raises(ValueError, match="looks like a FASTG file"):
        gfa2fastg(str(fastg_file), str(tmp_path / "converted_graph.fastg"))


def test_gfa2fastg_logs_error_for_fastg_input(tmp_path, monkeypatch):
    fastg_file = tmp_path / "graph.fastg"
    fastg_file.write_text(">A:B';\nACGT\n>B';\nTTTT\n")
    logged = []

    monkeypatch.setattr(
        gfa2fastg_module.logger, "error", lambda message: logged.append(message)
    )

    with pytest.raises(ValueError, match="looks like a FASTG file"):
        gfa2fastg(str(fastg_file), str(tmp_path / "converted_graph.fastg"))

    assert logged == [
        "graph.fastg looks like a FASTG file. "
        "The gfa2fastg subcommand expects GFA input."
    ]
