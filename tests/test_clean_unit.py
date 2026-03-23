#!/usr/bin/env python3

import importlib

from agtools.commands.clean import _write_filtered_graph, clean

clean_module = importlib.import_module("agtools.commands.clean")


class DummyParser:
    def __init__(self, index, sequences):
        self.index = index
        self._sequences = sequences

    def get_sequence(self, seq_name):
        return self._sequences[seq_name]


def test_write_filtered_graph_fills_missing_sequences_and_prunes_edges(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "S\ts1\t\tLN:i:4\n"
        "S\ts2\tTT\n"
        "L\ts1\t+\ts2\t-\t1M\n"
        "P\tp1\ts1+,s2-\t*\n"
        "P\tp2\ts1+\t*\n"
        "W\tw1\t*\t*\t*\t>s1<s2\n"
        "W\tw2\t*\t*\t*\t>s1\n"
    )

    parser = DummyParser(index={"s1": 0}, sequences={"s1": "ACGT"})

    output_file = _write_filtered_graph({"s2"}, parser, str(gfa_file), str(tmp_path))
    content = (tmp_path / "cleaned_graph.gfa").read_text()

    assert output_file == str(tmp_path / "cleaned_graph.gfa")
    assert "S\ts1\tACGT\tLN:i:4" in content
    assert "S\ts2\tTT" not in content
    assert "L\ts1\t+\ts2\t-\t1M" not in content
    assert "P\tp1\ts1+,s2-\t*" not in content
    assert "W\tw1\t*\t*\t*\t>s1<s2" not in content
    assert "P\tp2\ts1+\t*" in content
    assert "W\tw2\t*\t*\t*\t>s1" in content


def test_clean_removes_segments_missing_from_fasta_index(tmp_path, monkeypatch):
    class DummyUG:
        segment_names = ["s1", "s2"]

    parser = DummyParser(index={"s1": 0}, sequences={"s1": "ACGT"})

    monkeypatch.setattr(
        clean_module.UnitigGraph, "from_gfa", staticmethod(lambda _: DummyUG())
    )
    monkeypatch.setattr(clean_module, "FastaParser", lambda *_args, **_kwargs: parser)

    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\ts1\t\tLN:i:4\nS\ts2\tTT\nP\tp1\ts1+,s2-\t*\n")

    output_file = clean(str(gfa_file), fasta="contigs.fasta", assembler="general", output_path=str(tmp_path))
    content = (tmp_path / "cleaned_graph.gfa").read_text()

    assert output_file == str(tmp_path / "cleaned_graph.gfa")
    assert "S\ts1\tACGT\tLN:i:4" in content
    assert "S\ts2\tTT" not in content
    assert "P\tp1\ts1+,s2-\t*" not in content
