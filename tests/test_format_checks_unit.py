#!/usr/bin/env python3

import importlib

import pytest

from agtools.commands._format_checks import (
    _scan_graph_file,
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


def test_scan_graph_file_skips_blank_and_comment_lines_and_tracks_unknown_first_tag(
    tmp_path,
):
    graph_file = tmp_path / "graph.txt"
    graph_file.write_text("\n# comment\nZZ\tmystery\nVT\tseg1\tAAAA\n")

    first_format, counts, file_name = _scan_graph_file(str(graph_file))

    assert first_format == "ZZ"
    assert counts == {"gfa_segments": 0, "fastg_headers": 0, "asqg_segments": 1}
    assert file_name == "graph.txt"


def test_validate_gfa_input_rejects_asqg_file(tmp_path):
    graph_file = tmp_path / "graph.asqg"
    graph_file.write_text("VT\tseg1\tAAAA\n")

    with pytest.raises(ValueError, match="looks like an ASQG file"):
        validate_gfa_input(str(graph_file), "gfa2fasta")


def test_validate_gfa_input_rejects_unknown_non_gfa_file(tmp_path):
    graph_file = tmp_path / "graph.txt"
    graph_file.write_text("\n# comment\nZZ\tmystery\n")

    with pytest.raises(ValueError, match="No GFA segments were found"):
        validate_gfa_input(str(graph_file), "gfa2fasta")


def test_validate_fastg_input_rejects_asqg_file(tmp_path):
    graph_file = tmp_path / "graph.asqg"
    graph_file.write_text("VT\tseg1\tAAAA\n")

    with pytest.raises(ValueError, match="looks like an ASQG file"):
        validate_fastg_input(str(graph_file), "fastg2gfa")


def test_validate_fastg_input_rejects_unknown_non_fastg_file(tmp_path):
    graph_file = tmp_path / "graph.txt"
    graph_file.write_text("\n# comment\nZZ\tmystery\n")

    with pytest.raises(ValueError, match="No FASTG segment headers were found"):
        validate_fastg_input(str(graph_file), "fastg2gfa")


def test_validate_asqg_input_rejects_unknown_non_asqg_file(tmp_path):
    graph_file = tmp_path / "graph.txt"
    graph_file.write_text("\n# comment\nZZ\tmystery\n")

    with pytest.raises(ValueError, match="No ASQG segments were found"):
        validate_asqg_input(str(graph_file), "asqg2gfa")
