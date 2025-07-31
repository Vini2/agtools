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


def test_concat_seg_error(runner, tmp_dir):
    outpath = tmp_dir
    graph_1 = DATADIR / "test_graph.gfa"
    graph_2 = DATADIR / "test_graph.gfa"
    args = f"-g {graph_1} -g {graph_2} -o {outpath}".split()
    r = runner.invoke(concat, args, catch_exceptions=False)
    print(r.exit_code)
    assert r.exit_code == 1, r.output


def test_concat_path_error(runner, tmp_dir):
    outpath = tmp_dir
    graph_1 = DATADIR / "test_graph.gfa"
    graph_2 = DATADIR / "test_graph_2.gfa"
    args = f"-g {graph_1} -g {graph_2} -o {outpath}".split()
    r = runner.invoke(concat, args, catch_exceptions=False)
    print(r.exit_code)
    assert r.exit_code == 1, r.output


def test_concat_walk_error(runner, tmp_dir):
    outpath = tmp_dir
    graph_1 = DATADIR / "test_graph.gfa"
    graph_2 = DATADIR / "test_graph_3.gfa"
    args = f"-g {graph_1} -g {graph_2} -o {outpath}".split()
    r = runner.invoke(concat, args, catch_exceptions=False)
    print(r.exit_code)
    assert r.exit_code == 1, r.output
