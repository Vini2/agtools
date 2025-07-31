import subprocess

import pytest

__author__ = "Vijini Mallawaarachchi"
__credits__ = ["Vijini Mallawaarachchi"]


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    """set the working directory for all tests"""
    monkeypatch.chdir(tmp_dir)


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0

    Parameters
    ----------

    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_agtools():
    """test agtools"""
    cmd = "agtools --help"
    exec_command(cmd)


def test_agtools_version():
    """test agtools version"""
    cmd = "agtools --version"
    exec_command(cmd)


def test_agtools_stats():
    """test agtools stats"""
    cmd = "agtools stats --help"
    exec_command(cmd)


def test_agtools_rename():
    """test agtools rename"""
    cmd = "agtools rename --help"
    exec_command(cmd)


def test_agtools_concatenate():
    """test agtools concatenate"""
    cmd = "agtools concat --help"
    exec_command(cmd)


def test_agtools_filter():
    """test agtools filter"""
    cmd = "agtools filter --help"
    exec_command(cmd)


def test_agtools_clean():
    """test agtools clean"""
    cmd = "agtools clean --help"
    exec_command(cmd)


def test_agtools_component():
    """test agtools component"""
    cmd = "agtools component --help"
    exec_command(cmd)


def test_agtools_fastg2gfa():
    """test agtools fastg2gfa"""
    cmd = "agtools fastg2gfa --help"
    exec_command(cmd)


def test_agtools_asqg2gfa():
    """test agtools asqg2gfa"""
    cmd = "agtools asqg2gfa --help"
    exec_command(cmd)


def test_agtools_gfa2fastg():
    """test agtools gfa2fastg"""
    cmd = "agtools gfa2fastg --help"
    exec_command(cmd)


def test_agtools_gfa2dot():
    """test agtools gfa2dot"""
    cmd = "agtools gfa2dot --help"
    exec_command(cmd)


def test_agtools_gfa2fasta():
    """test agtools gfa2fasta"""
    cmd = "agtools gfa2fasta --help"
    exec_command(cmd)


def test_agtools_gfa2adj():
    """test agtools gfa2adj"""
    cmd = "agtools gfa2adj --help"
    exec_command(cmd)
