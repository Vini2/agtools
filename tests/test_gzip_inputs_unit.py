#!/usr/bin/env python3
"""End-to-end checks that gzip-compressed graph inputs are accepted."""

import gzip

import pytest

from agtools.commands.asqg2gfa import asqg2gfa
from agtools.commands.fastg2gfa import fastg2gfa
from agtools.commands.gfa2asqg import gfa2asqg
from agtools.commands.gfa2fasta import gfa2fasta
from agtools.commands.gfa2fastg import gfa2fastg
from agtools.core.unitig_graph import UnitigGraph

GFA = "H\tVN:Z:1.0\nS\tseg1\tATGCATGCAT\nS\tseg2\tGGTTCCAAGG\nL\tseg1\t+\tseg2\t+\t0M\n"
FASTG = ">seg1:seg2;\nATGCATGCAT\n>seg2;\nGGTTCCAAGG\n"
ASQG = (
    "HT\tVN:i:1\n"
    "VT\tseg1\tATGCATGCAT\n"
    "VT\tseg2\tGGTTCCAAGG\n"
    "ED\tseg1 seg2 0 9 10 0 9 10 0 0\n"
)


def _gzip_write(path, text):
    with gzip.open(path, "wt") as handle:
        handle.write(text)
    return path


@pytest.fixture
def gfa_pair(tmp_path):
    """A GFA file and its gzip-compressed twin."""
    plain = tmp_path / "graph.gfa"
    plain.write_text(GFA)
    return plain, _gzip_write(tmp_path / "graph.gfa.gz", GFA)


def test_unitig_graph_loads_from_gzipped_gfa(gfa_pair):
    plain, compressed = gfa_pair

    from_plain = UnitigGraph.from_gfa(str(plain))
    from_gzip = UnitigGraph.from_gfa(str(compressed))

    assert from_gzip.vcount == from_plain.vcount == 2
    assert from_gzip.lcount == from_plain.lcount == 1
    assert from_gzip.segment_names == from_plain.segment_names


def test_segment_sequence_lookup_works_on_gzipped_gfa(gfa_pair):
    _, compressed = gfa_pair

    graph = UnitigGraph.from_gfa(str(compressed))
    sequence = graph.get_segment_sequence("seg2")

    assert sequence == "GGTTCCAAGG"
    assert isinstance(sequence, str)


def test_gfa2fasta_accepts_gzipped_input(gfa_pair, tmp_path):
    plain, compressed = gfa_pair

    from_plain = tmp_path / "plain.fasta"
    from_gzip = tmp_path / "gzip.fasta"
    gfa2fasta(str(plain), str(from_plain))
    gfa2fasta(str(compressed), str(from_gzip))

    assert from_gzip.read_text() == from_plain.read_text()
    assert ">seg1\nATGCATGCAT\n" in from_gzip.read_text()


def test_gfa2fastg_accepts_gzipped_input(gfa_pair, tmp_path):
    plain, compressed = gfa_pair

    from_plain = tmp_path / "plain.fastg"
    from_gzip = tmp_path / "gzip.fastg"
    gfa2fastg(str(plain), str(from_plain))
    gfa2fastg(str(compressed), str(from_gzip))

    assert from_gzip.read_text() == from_plain.read_text()


def test_gfa2asqg_accepts_gzipped_input(gfa_pair, tmp_path):
    plain, compressed = gfa_pair

    from_plain = tmp_path / "plain.asqg"
    from_gzip = tmp_path / "gzip.asqg"
    gfa2asqg(str(plain), str(from_plain))
    gfa2asqg(str(compressed), str(from_gzip))

    assert from_gzip.read_text() == from_plain.read_text()


def test_fastg2gfa_accepts_gzipped_input(tmp_path):
    plain = tmp_path / "graph.fastg"
    plain.write_text(FASTG)
    compressed = _gzip_write(tmp_path / "graph.fastg.gz", FASTG)

    from_plain = tmp_path / "plain.gfa"
    from_gzip = tmp_path / "gzip.gfa"
    fastg2gfa(str(plain), 10, str(from_plain))
    fastg2gfa(str(compressed), 10, str(from_gzip))

    assert from_gzip.read_text() == from_plain.read_text()


def test_asqg2gfa_accepts_gzipped_input(tmp_path):
    plain = tmp_path / "graph.asqg"
    plain.write_text(ASQG)
    compressed = _gzip_write(tmp_path / "graph.asqg.gz", ASQG)

    from_plain = tmp_path / "plain.gfa"
    from_gzip = tmp_path / "gzip.gfa"
    asqg2gfa(str(plain), str(from_plain))
    asqg2gfa(str(compressed), str(from_gzip))

    assert from_gzip.read_text() == from_plain.read_text()


def test_format_validation_reads_gzipped_files(tmp_path):
    compressed = _gzip_write(tmp_path / "graph.fastg.gz", FASTG)

    with pytest.raises(ValueError, match="looks like a FASTG file"):
        gfa2fasta(str(compressed), str(tmp_path / "out.fasta"))
