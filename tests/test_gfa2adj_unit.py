#!/usr/bin/env python3

import importlib

import pandas as pd
import pytest

from agtools.commands.gfa2adj import gfa2adj

gfa2adj_module = importlib.import_module("agtools.commands.gfa2adj")


def test_gfa2adj_writes_comma_delimited_csv(tmp_path, monkeypatch):
    class DummyUG:
        def get_adjacency_matrix(self, type="matrix"):
            assert type == "pandas"
            return pd.DataFrame([[0, 1], [1, 0]], index=["s1", "s2"], columns=["s1", "s2"])

    monkeypatch.setattr(
        gfa2adj_module.UnitigGraph, "from_gfa", staticmethod(lambda _path: DummyUG())
    )
    monkeypatch.setattr(gfa2adj_module, "validate_gfa_input", lambda *_args: None)

    target = tmp_path / "adjacency_matrix.csv"
    output_file = gfa2adj("graph.gfa", delimiter="comma", output_path=str(target))
    content = target.read_text().splitlines()[0]

    assert output_file == str(tmp_path / "adjacency_matrix.csv")
    assert content == ",s1,s2"


def test_gfa2adj_writes_tab_delimited_tsv(tmp_path, monkeypatch):
    class DummyUG:
        def get_adjacency_matrix(self, type="matrix"):
            assert type == "pandas"
            return pd.DataFrame([[0, 1], [1, 0]], index=["s1", "s2"], columns=["s1", "s2"])

    monkeypatch.setattr(
        gfa2adj_module.UnitigGraph, "from_gfa", staticmethod(lambda _path: DummyUG())
    )
    monkeypatch.setattr(gfa2adj_module, "validate_gfa_input", lambda *_args: None)

    target = tmp_path / "adjacency_matrix.tsv"
    output_file = gfa2adj("graph.gfa", delimiter="tab", output_path=str(target))
    content = target.read_text().splitlines()[0]

    assert output_file == str(tmp_path / "adjacency_matrix.tsv")
    assert content == "\ts1\ts2"


def test_gfa2adj_rejects_fastg_input(tmp_path):
    fastg_file = tmp_path / "graph.fastg"
    fastg_file.write_text(">A:B';\nACGT\n>B';\nTTTT\n")

    with pytest.raises(ValueError, match="looks like a FASTG file"):
        gfa2adj(str(fastg_file), delimiter="comma", output_path=str(tmp_path / "out.csv"))
