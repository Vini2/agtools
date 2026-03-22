#!/usr/bin/env python3

import gzip

import pytest

from agtools.core.fasta_parser import FastaParser


def test_plain_fasta_index_and_sequence_lookup(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1 description\nATGC\nTTAA\n>seq2\nGGCC\n")

    parser = FastaParser(fasta_file)

    assert set(parser.index.keys()) == {"seq1", "seq2"}
    assert str(parser.get_sequence("seq1")) == "ATGCTTAA"
    assert parser.get_index("seq2") == parser.index["seq2"]


def test_missing_sequence_warns_and_returns_empty_string(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1\nATGC\n")

    parser = FastaParser(fasta_file)

    with pytest.warns(RuntimeWarning, match="not found"):
        sequence = parser.get_sequence("unknown_seq")

    assert sequence == ""


def test_get_index_raises_for_unknown_sequence(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1\nATGC\n")

    parser = FastaParser(fasta_file)

    with pytest.raises(KeyError, match="unknown_seq not found in the index"):
        parser.get_index("unknown_seq")


def test_myloasm_header_parsing_uses_prefix_before_underscore(tmp_path):
    fasta_file = tmp_path / "myloasm.fasta"
    fasta_file.write_text(">u913838ctg_42\nATGCATGC\n")

    parser = FastaParser(fasta_file, assembler="myloasm")

    assert "u913838ctg" in parser.index
    assert str(parser.get_sequence("u913838ctg")) == "ATGCATGC"


def test_megahit_mapping_is_used_for_sequence_lookup(tmp_path):
    fasta_file = tmp_path / "megahit.fasta"
    fasta_file.write_text(">k141_4704\nAATTCCGG\n")

    parser = FastaParser(
        fasta_file, assembler="megahit", mapping={"NODE_1_length_205_cov_1.0000_ID_1": "k141_4704"}
    )

    assert (
        str(parser.get_sequence("NODE_1_length_205_cov_1.0000_ID_1")) == "AATTCCGG"
    )


def test_gzipped_fasta_is_supported(tmp_path):
    fasta_file = tmp_path / "contigs.fasta.gz"
    with gzip.open(fasta_file, "wt") as f:
        f.write(">seq1\nATGC\n>seq2\nGGTTAA\n")

    parser = FastaParser(fasta_file)

    assert parser.gzipped is True
    assert str(parser.get_sequence("seq2")) == "GGTTAA"
