#!/usr/bin/env python3

import pathlib

import pytest
from click.testing import CliRunner

from agtools.cli import *

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


def _assert_output_file_created(path):
    output_path = pathlib.Path(str(path))
    assert output_path.exists()
    assert output_path.is_file()
    assert output_path.stat().st_size > 0


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    """set the working directory for all tests"""
    monkeypatch.chdir(tmp_dir)


@pytest.fixture(scope="session")
def runner():
    """exportrc works correctly."""
    return CliRunner()


def test_agtools_stats(runner, tmp_dir):
    outpath = tmp_dir / "stats" / "graph_stats.txt"
    graph = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(stats, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_rename_seg(runner, tmp_dir):
    outpath = tmp_dir / "rename_seg" / "renamed_graph.gfa"
    graph = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    prefix = "test"
    args = f"-g {graph} -p {prefix} -o {outpath}".split()
    r = runner.invoke(rename, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_rename_path(runner, tmp_dir):
    outpath = tmp_dir / "rename_path" / "renamed_graph.gfa"
    graph = DATADIR / "test_path.gfa"
    prefix = "test"
    args = f"-g {graph} -p {prefix} -o {outpath}".split()
    r = runner.invoke(rename, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_rename_walk(runner, tmp_dir):
    outpath = tmp_dir / "rename_walk" / "renamed_graph.gfa"
    graph = DATADIR / "test_walk.gfa"
    prefix = "test"
    args = f"-g {graph} -p {prefix} -o {outpath}".split()
    r = runner.invoke(rename, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_concat(runner, tmp_dir):
    outpath = tmp_dir / "concat" / "concatenated_graph.gfa"
    graph_1 = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    graph_2 = DATADIR / "test_graph.gfa"
    args = f"-g {graph_1} -g {graph_2} -o {outpath}".split()
    r = runner.invoke(concat, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_filter(runner, tmp_dir):
    outpath = tmp_dir / "filter" / "filtered_graph.gfa"
    graph = DATADIR / "test_graph.gfa"
    min_length = 1000
    args = f"-g {graph} -l {min_length} -o {outpath}".split()
    r = runner.invoke(filter, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_clean(runner, tmp_dir):
    outpath = tmp_dir / "clean" / "cleaned_graph.gfa"
    graph = DATADIR / "test_graph_1.gfa"
    fasta = DATADIR / "test_fasta.fasta"
    args = f"-g {graph} -f {fasta} -a general -o {outpath}".split()
    r = runner.invoke(clean, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_clean_myloasm(runner, tmp_dir):
    outpath = tmp_dir / "clean_myloasm" / "cleaned_graph.gfa"
    graph = DATADIR / "myloasm" / "final_contig_graph.gfa"
    fasta = DATADIR / "myloasm" / "assembly_primary.fa"
    args = f"-g {graph} -f {fasta} -a myloasm -o {outpath}".split()
    r = runner.invoke(clean, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_component(runner, tmp_dir):
    outpath = tmp_dir / "component" / "component_graph.gfa"
    graph = DATADIR / "test_graph.gfa"
    segment = "seg4"
    args = f"-g {graph} -s {segment} -o {outpath}".split()
    r = runner.invoke(component, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_fastg2gfa(runner, tmp_dir):
    outpath = tmp_dir / "fastg2gfa" / "converted_graph.gfa"
    graph = DATADIR / "final.graph.fastg"
    k = 141
    args = f"-g {graph} -k {k} -o {outpath}".split()
    r = runner.invoke(fastg2gfa, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_fastg2gfa_rejects_gfa_input(runner, tmp_dir):
    outpath = tmp_dir / "fastg2gfa_error" / "converted_graph.gfa"
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -k 141 -o {outpath}".split()
    r = runner.invoke(fastg2gfa, args, catch_exceptions=False)
    assert r.exit_code == 1
    assert "looks like a GFA file" in r.output


def test_agtools_gfa2fastg(runner, tmp_dir):
    outpath = tmp_dir / "gfa2fastg" / "converted_graph.fastg"
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2fastg, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_gfa2fastg_rejects_fastg_input(runner, tmp_dir):
    outpath = tmp_dir / "gfa2fastg_error" / "converted_graph.fastg"
    graph = DATADIR / "final.graph.fastg"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2fastg, args, catch_exceptions=False)
    assert r.exit_code == 1
    assert "looks like a FASTG file" in r.output


def test_agtools_asqg2gfa(runner, tmp_dir):
    outpath = tmp_dir / "asqg2gfa" / "converted_graph.gfa"
    graph = DATADIR / "ESC" / "default-graph.asqg"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(asqg2gfa, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_asqg2gfa_rejects_gfa_input(runner, tmp_dir):
    outpath = tmp_dir / "asqg2gfa_error" / "converted_graph.gfa"
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(asqg2gfa, args, catch_exceptions=False)
    assert r.exit_code == 1
    assert "looks like a GFA file" in r.output


def test_agtools_gfa2asqg(runner, tmp_dir):
    outpath = tmp_dir / "gfa2asqg" / "converted_graph.asqg"
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2asqg, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_gfa2asqg_rejects_asqg_input(runner, tmp_dir):
    outpath = tmp_dir / "gfa2asqg_error" / "converted_graph.asqg"
    graph = DATADIR / "ESC" / "default-graph.asqg"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2asqg, args, catch_exceptions=False)
    assert r.exit_code == 1
    assert "looks like an ASQG file" in r.output


def test_agtools_gfa2dot(runner, tmp_dir):
    outpath = tmp_dir / "gfa2dot" / "graph.dot"
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2dot, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_gfa2dot_rejects_fastg_input(runner, tmp_dir):
    outpath = tmp_dir / "gfa2dot_error" / "graph.dot"
    graph = DATADIR / "final.graph.fastg"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2dot, args, catch_exceptions=False)
    assert r.exit_code == 1
    assert "looks like a FASTG file" in r.output


def test_agtools_gfa2dot_abyss(runner, tmp_dir):
    outpath = tmp_dir / "gfa2dot_abyss" / "graph.gv"
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -ab -o {outpath}".split()
    r = runner.invoke(gfa2dot, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_gfa2fasta(runner, tmp_dir):
    outpath = tmp_dir / "gfa2fasta" / "segments.fasta"
    graph = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2fasta, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_gfa2fasta_rejects_fastg_input(runner, tmp_dir):
    outpath = tmp_dir / "gfa2fasta_error" / "segments.fasta"
    graph = DATADIR / "final.graph.fastg"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2fasta, args, catch_exceptions=False)
    assert r.exit_code == 1
    assert "looks like a FASTG file" in r.output


def test_agtools_gfa2adj(runner, tmp_dir):
    outpath = tmp_dir / "gfa2adj" / "adjacency_matrix.csv"
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2adj, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    _assert_output_file_created(outpath)


def test_agtools_gfa2adj_rejects_fastg_input(runner, tmp_dir):
    outpath = tmp_dir / "gfa2adj_error" / "adjacency_matrix.csv"
    graph = DATADIR / "final.graph.fastg"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2adj, args, catch_exceptions=False)
    assert r.exit_code == 1
    assert "looks like a FASTG file" in r.output


def test_agtools_no_log_file_by_default(runner, tmp_dir):
    outpath = tmp_dir / "no_log" / "adjacency_matrix.csv"
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2adj, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    assert not (tmp_dir / "agtools.log").exists()


def test_agtools_log_file_is_optional(runner, tmp_dir):
    outpath = tmp_dir / "with_log" / "adjacency_matrix.csv"
    graph = DATADIR / "test_graph.gfa"
    log_file = tmp_dir / "custom-agtools.log"
    args = f"-g {graph} -o {outpath} --log-file {log_file}".split()
    r = runner.invoke(gfa2adj, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    assert log_file.exists()
