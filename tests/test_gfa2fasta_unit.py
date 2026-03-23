#!/usr/bin/env python3

from agtools.commands.gfa2fasta import (
    _get_segment_sequences,
    _write_segment_sequences,
    gfa2fasta,
)


def test_get_segment_sequences_normalizes_non_nucleotide_chars(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "H\tVN:Z:1.0\n" "S\tseg1\tatgcnN-*x\n" "S\tseg2\tGGttac\n" "L\tseg1\t+\tseg2\t+\t1M\n"
    )

    sequences = _get_segment_sequences(str(gfa_file))

    assert [record.id for record in sequences] == ["seg1", "seg2"]
    assert str(sequences[0].seq) == "ATGC"
    assert str(sequences[1].seq) == "GGTTAC"


def test_write_segment_sequences_outputs_fasta(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tATGC\n")
    sequences = _get_segment_sequences(str(gfa_file))

    output_file = _write_segment_sequences(sequences, str(tmp_path))
    content = (tmp_path / "segments.fasta").read_text()

    assert output_file == str(tmp_path / "segments.fasta")
    assert ">seg1" in content
    assert "ATGC" in content


def test_gfa2fasta_end_to_end(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tseg1\tATGC\nS\tseg2\tGGTT\n")

    output_file = gfa2fasta(str(gfa_file), str(tmp_path))
    content = (tmp_path / "segments.fasta").read_text()

    assert output_file == str(tmp_path / "segments.fasta")
    assert ">seg1" in content
    assert ">seg2" in content

