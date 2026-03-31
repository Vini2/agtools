#!/usr/bin/env python3

import importlib

import pytest

from agtools.commands._format_checks import (
    validate_asqg_input,
    validate_fastg_input,
    validate_gfa_input,
)

format_checks_module = importlib.import_module("agtools.commands._format_checks")


def test_validate_gfa_input_logs_and_raises_for_fastg_file(tmp_path, monkeypatch):
    graph_file = tmp_path / "graph.fastg"
    graph_file.write_text(">A:B';\nACGT\n>B';\nTTTT\n")
    logged = []

    monkeypatch.setattr(
        format_checks_module.logger, "error", lambda message: logged.append(message)
    )

    with pytest.raises(ValueError, match="looks like a FASTG file"):
        validate_gfa_input(str(graph_file), "gfa2fastg")

    assert logged == [
        "graph.fastg looks like a FASTG file. "
        "The gfa2fastg subcommand expects GFA input."
    ]


def test_validate_fastg_input_rejects_gfa_file(tmp_path):
    graph_file = tmp_path / "graph.gfa"
    graph_file.write_text("S\tseg1\tATGC\n")

    with pytest.raises(ValueError, match="looks like a GFA file"):
        validate_fastg_input(str(graph_file), "fastg2gfa")


def test_validate_asqg_input_rejects_fastg_file(tmp_path):
    graph_file = tmp_path / "graph.fastg"
    graph_file.write_text(">A:B';\nACGT\n>B';\nTTTT\n")

    with pytest.raises(ValueError, match="looks like a FASTG file"):
        validate_asqg_input(str(graph_file), "asqg2gfa")
