#!/usr/bin/env python3

import gzip

import pytest

from agtools.core._io import (
    is_bgzf,
    is_gzip,
    open_text,
    strip_compression_suffix,
)


def test_is_gzip_detects_by_magic_bytes(tmp_path):
    plain = tmp_path / "graph.gfa"
    plain.write_text("S\tseg1\tATGC\n")

    compressed = tmp_path / "graph.gfa.gz"
    with gzip.open(compressed, "wt") as handle:
        handle.write("S\tseg1\tATGC\n")

    # A gzip file without a .gz suffix is still detected.
    disguised_gz = tmp_path / "nosuffix"
    disguised_gz.write_bytes(compressed.read_bytes())

    assert is_gzip(plain) is False
    assert is_gzip(compressed) is True
    assert is_gzip(disguised_gz) is True
    assert is_gzip(tmp_path / "missing.gfa") is False


def test_is_bgzf_distinguishes_bgzf_from_plain_gzip(tmp_path, write_bgzf):
    plain = tmp_path / "a.gfa"
    plain.write_text("S\tseg1\tATGC\n")

    ordinary = tmp_path / "b.gfa.gz"
    with gzip.open(ordinary, "wt") as handle:
        handle.write("S\tseg1\tATGC\n")

    blocked = write_bgzf(tmp_path / "c.gfa.gz", "S\tseg1\tATGC\n")

    assert is_bgzf(plain) is False
    assert is_bgzf(ordinary) is False
    assert is_bgzf(blocked) is True


def test_open_text_reads_plain_and_gzip(tmp_path):
    plain = tmp_path / "a.gfa"
    plain.write_text("line1\nline2\n")

    compressed = tmp_path / "b.gfa.gz"
    with gzip.open(compressed, "wt") as handle:
        handle.write("line1\nline2\n")

    with open_text(plain) as handle:
        assert handle.read() == "line1\nline2\n"

    with open_text(compressed) as handle:
        assert handle.read() == "line1\nline2\n"


def test_open_text_supports_seeking_in_gzip(tmp_path):
    compressed = tmp_path / "b.gfa.gz"
    with gzip.open(compressed, "wt") as handle:
        handle.write("line1\nline2\n")

    with open_text(compressed) as handle:
        handle.readline()
        offset = handle.tell()
        handle.readline()
        handle.seek(offset)
        assert handle.readline() == "line2\n"


def test_open_text_compresses_on_write_when_path_ends_in_gz(tmp_path):
    target = tmp_path / "out.txt.gz"

    with open_text(target, "wt") as handle:
        handle.write("hello\n")

    assert is_gzip(target) is True
    with gzip.open(target, "rt") as handle:
        assert handle.read() == "hello\n"


def test_open_text_rejects_binary_modes(tmp_path):
    with pytest.raises(ValueError, match="text modes"):
        open_text(tmp_path / "x.gfa", "rb")


def test_strip_compression_suffix():
    assert strip_compression_suffix("contigs.fasta.gz") == "contigs.fasta"
    assert strip_compression_suffix("contigs.fasta") == "contigs.fasta"
