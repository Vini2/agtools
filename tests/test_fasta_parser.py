#!/usr/bin/env python3

import gzip
import os

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
        fasta_file,
        assembler="megahit",
        mapping={"NODE_1_length_205_cov_1.0000_ID_1": "k141_4704"},
    )

    assert str(parser.get_sequence("NODE_1_length_205_cov_1.0000_ID_1")) == "AATTCCGG"


def test_gzipped_fasta_is_supported(tmp_path):
    fasta_file = tmp_path / "contigs.fasta.gz"
    with gzip.open(fasta_file, "wt") as f:
        f.write(">seq1\nATGC\n>seq2\nGGTTAA\n")

    parser = FastaParser(fasta_file)

    assert parser.gzipped is True
    assert str(parser.get_sequence("seq2")) == "GGTTAA"


def test_plain_fasta_uses_the_rsbio_seq_index(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1\nATGC\n>seq2\nGGCC\n")

    parser = FastaParser(fasta_file)

    assert parser.indexed is True
    assert (tmp_path / "contigs.fasta.fai").exists()


def test_existing_index_is_reused_and_stale_index_is_rebuilt(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1\nATGC\n")
    FastaParser(fasta_file)

    fai_file = tmp_path / "contigs.fasta.fai"
    first_mtime = fai_file.stat().st_mtime

    FastaParser(fasta_file)
    assert fai_file.stat().st_mtime == first_mtime

    # A FASTA file newer than its index forces a rebuild.
    fasta_file.write_text(">seq1\nATGC\n>seq2\nTTTT\n")
    os.utime(fai_file, (first_mtime - 100, first_mtime - 100))

    parser = FastaParser(fasta_file)
    assert parser.get_sequence("seq2") == "TTTT"


def test_bgzf_fasta_is_indexed_for_random_access(tmp_path, write_bgzf):
    fasta_file = write_bgzf(tmp_path / "contigs.fa.gz", ">seq1\nATGC\n>seq2\nGGTTAA\n")

    parser = FastaParser(fasta_file)

    assert parser.gzipped is True
    assert parser.bgzf is True
    assert parser.indexed is True
    assert parser.get_sequence("seq2") == "GGTTAA"


def test_ordinary_gzip_falls_back_to_scanning(tmp_path):
    fasta_file = tmp_path / "contigs.fasta.gz"
    with gzip.open(fasta_file, "wt") as f:
        f.write(">seq1\nATGC\n>seq2\nGGTTAA\n")

    parser = FastaParser(fasta_file)

    assert parser.bgzf is False
    assert parser.indexed is False
    assert not (tmp_path / "contigs.fasta.gz.fai").exists()
    assert parser.get_sequence("seq1") == "ATGC"


def test_unrecognised_extension_falls_back_to_scanning(tmp_path):
    fasta_file = tmp_path / "contigs.txt"
    fasta_file.write_text(">seq1\nATGC\n")

    parser = FastaParser(fasta_file)

    assert parser.indexed is False
    assert parser.get_sequence("seq1") == "ATGC"


def test_sequences_are_returned_as_plain_strings(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1\nATGC\n")

    sequence = FastaParser(fasta_file).get_sequence("seq1")

    assert isinstance(sequence, str)


def test_iter_sequences_streams_every_record(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1 desc\nATGC\nTTAA\n>seq2\nGGCC\n")

    parser = FastaParser(fasta_file)

    assert list(parser.iter_sequences()) == [("seq1", "ATGCTTAA"), ("seq2", "GGCC")]


def test_iter_sequences_works_on_gzipped_input(tmp_path):
    fasta_file = tmp_path / "contigs.fasta.gz"
    with gzip.open(fasta_file, "wt") as f:
        f.write(">seq1\nATGC\n>seq2\nGGTTAA\n")

    parser = FastaParser(fasta_file)

    assert dict(parser.iter_sequences()) == {"seq1": "ATGC", "seq2": "GGTTAA"}


def test_iter_lengths_matches_sequence_lengths(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1 desc\nATGC\nTTAA\n>seq2\nGGCC\n")

    parser = FastaParser(fasta_file)

    assert parser.has_cheap_lengths() is True
    assert list(parser.iter_lengths()) == [("seq1", 8), ("seq2", 4)]
    assert dict(parser.iter_lengths()) == {
        name: len(seq) for name, seq in parser.iter_sequences()
    }


def test_iter_lengths_falls_back_without_an_index(tmp_path):
    fasta_file = tmp_path / "contigs.fasta.gz"
    with gzip.open(fasta_file, "wt") as f:
        f.write(">seq1\nATGC\n>seq2\nGGTTAA\n")

    parser = FastaParser(fasta_file)

    assert parser.has_cheap_lengths() is False
    assert dict(parser.iter_lengths()) == {"seq1": 4, "seq2": 6}


def test_index_is_lazy_and_does_not_materialise_offsets(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1\nATGC\n>seq2\nGGCC\n")

    parser = FastaParser(fasta_file)

    # Iteration, counting and lengths all stream off the .fai and retain
    # nothing, so none of them build the offset table.
    assert list(parser.index) == ["seq1", "seq2"]
    assert len(parser.index) == 2
    assert list(parser.iter_lengths()) == [("seq1", 4), ("seq2", 4)]
    assert parser.index._offsets is None

    # Membership and subscripting materialise it, once.
    assert "seq1" in parser.index
    assert parser.index["seq2"] == 17
    assert parser.index._offsets is not None


def test_indexed_reader_is_not_opened_until_a_sequence_is_read(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1\nATGC\n>seq2\nGGCC\n")

    parser = FastaParser(fasta_file)

    # Loading an index on the Rust side costs real memory on large
    # assemblies, so nothing opens it until a base is actually needed.
    assert parser._reader is None
    len(parser.index)
    list(parser.iter_lengths())
    assert parser._reader is None

    assert parser.get_sequence("seq1") == "ATGC"
    assert parser._reader is not None


def test_index_behaves_like_a_mapping(tmp_path):
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">seq1\nATGC\n>seq2\nGGCC\n")

    index = FastaParser(fasta_file).index

    assert set(index.keys()) == {"seq1", "seq2"}
    assert dict(index) == {"seq1": index["seq1"], "seq2": index["seq2"]}
    assert index.get("nope") is None
    with pytest.raises(KeyError):
        index["nope"]


def test_myloasm_index_deduplicates_shared_prefixes(tmp_path):
    fasta_file = tmp_path / "assembly_primary.fa"
    fasta_file.write_text(">ctgA_1\nATGC\n>ctgA_2\nGGGG\n>ctgB_1\nTTTT\n")

    parser = FastaParser(fasta_file, assembler="myloasm")

    assert list(parser.index) == ["ctgA", "ctgB"]
    assert parser.has_cheap_lengths() is False
    assert "ctgB" in parser.index
    assert parser.get_sequence("ctgB") == "TTTT"


def test_plain_gzip_warns_that_it_cannot_be_indexed(tmp_path, caplog):
    fasta_file = tmp_path / "contigs.fasta.gz"
    with gzip.open(fasta_file, "wt") as f:
        f.write(">seq1\nATGC\n")

    from agtools.log_config import logger

    messages = []
    handler_id = logger.add(messages.append, level="WARNING")
    try:
        parser = FastaParser(fasta_file)
    finally:
        logger.remove(handler_id)

    assert parser.indexed is False
    assert any("cannot be indexed" in str(m) for m in messages)
    assert any("bgzip" in str(m) for m in messages)


def test_bgzf_and_plain_do_not_warn(tmp_path, write_bgzf):
    from agtools.log_config import logger

    messages = []
    handler_id = logger.add(messages.append, level="WARNING")
    try:
        plain = tmp_path / "plain.fasta"
        plain.write_text(">seq1\nATGC\n")
        FastaParser(plain)
        FastaParser(write_bgzf(tmp_path / "blocked.fa.gz", ">seq1\nATGC\n"))
    finally:
        logger.remove(handler_id)

    assert messages == []
