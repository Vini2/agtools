#!/usr/bin/env python3

import pathlib

import pytest
from click.testing import CliRunner

from agtools.cli import *

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]

DATADIR = pathlib.Path(__file__).parent / "data"


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
    outpath = tmp_dir
    graph = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(stats, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_rename_seg(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    prefix = "test"
    args = f"-g {graph} -p {prefix} -o {outpath}".split()
    r = runner.invoke(rename, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_rename_path(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "test_path.gfa"
    prefix = "test"
    args = f"-g {graph} -p {prefix} -o {outpath}".split()
    r = runner.invoke(rename, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_rename_walk(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "test_walk.gfa"
    prefix = "test"
    args = f"-g {graph} -p {prefix} -o {outpath}".split()
    r = runner.invoke(rename, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_concat(runner, tmp_dir):
    outpath = tmp_dir
    graph_1 = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    graph_2 = DATADIR / "test_graph.gfa"
    args = f"-g {graph_1} -g {graph_2} -o {outpath}".split()
    r = runner.invoke(concat, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_filter(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "test_graph.gfa"
    min_length = 1000
    args = f"-g {graph} -l {min_length} -o {outpath}".split()
    r = runner.invoke(filter, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_clean(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "myloasm" / "final_contig_graph.gfa"
    fasta = DATADIR / "myloasm" / "assembly_primary.fa"
    args = f"-g {graph} -f {fasta} -a myloasm -o {outpath}".split()
    r = runner.invoke(clean, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_component(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "test_graph.gfa"
    segment = "seg4"
    args = f"-g {graph} -s {segment} -o {outpath}".split()
    r = runner.invoke(component, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_fastg2gfa(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "final.graph.fastg"
    k = 141
    args = f"-g {graph} -k {k} -o {outpath}".split()
    r = runner.invoke(fastg2gfa, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_asqg2gfa(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "ESC" / "default-graph.asqg"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(asqg2gfa, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_gfa2fastg(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2fastg, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_gfa2dot(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2dot, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_gfa2dot_abyss(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -ab -o {outpath}".split()
    r = runner.invoke(gfa2dot, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_gfa2fasta(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "ESC" / "assembly_graph_with_scaffolds.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2fasta, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_agtools_gfa2adj(runner, tmp_dir):
    outpath = tmp_dir
    graph = DATADIR / "test_graph.gfa"
    args = f"-g {graph} -o {outpath}".split()
    r = runner.invoke(gfa2adj, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
