#!/usr/bin/env python3

from agtools.commands.asqg2gfa import asqg2gfa
from agtools.commands.clean import clean
from agtools.commands.component import component
from agtools.commands.concat import concat
from agtools.commands.fastg2gfa import fastg2gfa
from agtools.commands.filter import filter
from agtools.commands.gfa2adj import gfa2adj
from agtools.commands.gfa2dot import gfa2dot
from agtools.commands.gfa2fasta import gfa2fasta
from agtools.commands.rename import rename
from agtools.commands.stats import stats


def _write_base_gfa(path):
    path.write_text(
        "H\tVN:Z:1.0\n"
        "S\ts1\tAAAA\n"
        "S\ts2\tCC\n"
        "S\ts3\tAT\n"
        "L\ts1\t+\ts2\t+\t1M\n"
        "P\tp1\ts1+,s2+\t*\n"
        "W\tw1\t*\t*\t*\t>s1<s2\n"
    )


def test_stats_output_file_content(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    _write_base_gfa(gfa_file)

    target = tmp_path / "graph_stats.txt"
    output_file = stats(str(gfa_file), str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "graph_stats.txt")
    assert "Number of segments: 3" in content
    assert "Number of links: 1" in content
    assert "Number of connected components: 2" in content
    assert "Number of self-loops: 0" in content
    assert "Total length of segments: 8 bp" in content
    assert "N50: 4 bp" in content
    assert "L50: 1 segment(s)" in content
    assert "GC content: 25.00%" in content


def test_rename_output_file_content(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    _write_base_gfa(gfa_file)

    target = tmp_path / "renamed_graph.gfa"
    output_file = rename(str(gfa_file), "pref", str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "renamed_graph.gfa")
    assert "S\tpref_s1\tAAAA" in content
    assert "L\tpref_s1\t+\tpref_s2\t+\t1M" in content
    assert "P\tpref_p1\tpref_s1+,pref_s2+\t*" in content
    assert "W\tpref_w1\t*\t*\t*\t>pref_s1<pref_s2" in content


def test_concat_output_file_content(tmp_path):
    gfa_1 = tmp_path / "g1.gfa"
    gfa_2 = tmp_path / "g2.gfa"
    gfa_1.write_text("H\tVN:Z:1.0\nS\ts1\tAAAA\n")
    gfa_2.write_text("S\ts2\tCC\nL\ts1\t+\ts2\t+\t1M\n")

    target = tmp_path / "concatenated_graph.gfa"
    output_file = concat([str(gfa_1), str(gfa_2)], str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "concatenated_graph.gfa")
    assert "H\tVN:Z:1.0" in content
    assert "S\ts1\tAAAA" in content
    assert "S\ts2\tCC" in content
    assert "L\ts1\t+\ts2\t+\t1M" in content


def test_filter_output_file_content(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    _write_base_gfa(gfa_file)

    target = tmp_path / "filtered_graph.gfa"
    output_file = filter(str(gfa_file), min_length=3, output_path=str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "filtered_graph.gfa")
    assert "S\ts1\tAAAA" in content
    assert "S\ts2\tCC" not in content
    assert "S\ts3\tAT" not in content
    assert "L\ts1\t+\ts2\t+\t1M" not in content
    assert "P\tp1\ts1+,s2+\t*" not in content
    assert "W\tw1\t*\t*\t*\t>s1<s2" not in content


def test_clean_output_file_content(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "S\ts1\t\tLN:i:4\n"
        "S\ts2\tCC\n"
        "L\ts1\t+\ts2\t+\t1M\n"
        "P\tp1\ts1+,s2+\t*\n"
    )
    fasta_file = tmp_path / "contigs.fasta"
    fasta_file.write_text(">s1\nAAAA\n")

    target = tmp_path / "cleaned_graph.gfa"
    output_file = clean(
        str(gfa_file), str(fasta_file), assembler="general", output_path=str(target)
    )
    content = target.read_text()

    assert output_file == str(tmp_path / "cleaned_graph.gfa")
    assert "S\ts1\tAAAA\tLN:i:4" in content
    assert "S\ts2\tCC" not in content
    assert "L\ts1\t+\ts2\t+\t1M" not in content
    assert "P\tp1\ts1+,s2+\t*" not in content


def test_component_output_file_content(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "S\ts1\tAAAA\n"
        "S\ts2\tCC\n"
        "S\ts3\tAT\n"
        "L\ts1\t+\ts2\t+\t1M\n"
        "P\tp1\ts1+,s2+\t*\n"
        "P\tp2\ts3+\t*\n"
    )

    target = tmp_path / "component_graph.gfa"
    output_file = component(str(gfa_file), segment="s1", output_path=str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "component_graph.gfa")
    assert "S\ts1\tAAAA" in content
    assert "S\ts2\tCC" in content
    assert "S\ts3\tAT" not in content
    assert "P\tp1\ts1+,s2+\t*" in content
    assert "P\tp2\ts3+\t*" not in content


def test_fastg2gfa_output_file_content(tmp_path):
    fastg_file = tmp_path / "graph.fastg"
    fastg_file.write_text(">A:B';\nATGC\n>B';\nGG\n")

    target = tmp_path / "converted_graph.gfa"
    output_file = fastg2gfa(str(fastg_file), k_overlap=55, gfa_path=str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "converted_graph.gfa")
    assert "S\tA\tATGC" in content
    assert "L\tA\t+\tB\t-\t55M" in content


def test_asqg2gfa_output_file_content(tmp_path):
    asqg_file = tmp_path / "graph.asqg"
    asqg_file.write_text(
        "VT\tcontig1\tAAAA\n"
        "VT\tcontig2\tCCCC\n"
        "ED\tcontig1 contig2 0 3 0 0 3 0 1\n"
    )

    target = tmp_path / "converted_graph.gfa"
    output_file = asqg2gfa(str(asqg_file), str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "converted_graph.gfa")
    assert "S\tcontig1\tAAAA" in content
    assert "S\tcontig2\tCCCC" in content
    assert "L\tcontig1\t+\tcontig2\t-\t3M" in content


def test_gfa2dot_outputs_file_for_both_formats(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    _write_base_gfa(gfa_file)

    dot_target = tmp_path / "graph.dot"
    abyss_target = tmp_path / "graph.gv"
    dot_file = gfa2dot(str(gfa_file), abyss=False, output_path=str(dot_target))
    abyss_file = gfa2dot(str(gfa_file), abyss=True, output_path=str(abyss_target))

    assert dot_file == str(dot_target)
    assert abyss_file == str(abyss_target)
    assert dot_target.exists()
    assert "digraph g {" in abyss_target.read_text()


def test_gfa2fasta_output_file_content(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\ts1\tatgcn-\nS\ts2\tGGTT\n")

    target = tmp_path / "segments.fasta"
    output_file = gfa2fasta(str(gfa_file), str(target))
    content = target.read_text()

    assert output_file == str(tmp_path / "segments.fasta")
    assert ">s1" in content
    assert "ATGC" in content
    assert ">s2" in content
    assert "GGTT" in content


def test_gfa2adj_output_file_content(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    _write_base_gfa(gfa_file)

    target = tmp_path / "adjacency_matrix.csv"
    output_file = gfa2adj(str(gfa_file), delimiter="comma", output_path=str(target))
    rows = target.read_text().splitlines()

    assert output_file == str(tmp_path / "adjacency_matrix.csv")
    assert rows[0] == ",s1,s2,s3"
    assert any(row.startswith("s1,0,1,0") for row in rows[1:])
    assert any(row.startswith("s2,1,0,0") for row in rows[1:])
    assert any(row.startswith("s3,0,0,0") for row in rows[1:])
