#!/usr/bin/env python3

import io

from agtools.commands._output import open_output_file, prepare_output_file


def test_prepare_output_file_creates_parent_directories_and_expands_user(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("HOME", str(tmp_path))

    output_file = prepare_output_file("~/nested/outputs/result.txt")

    assert output_file == str(tmp_path / "nested" / "outputs" / "result.txt")
    assert (tmp_path / "nested" / "outputs").is_dir()


def test_prepare_output_file_preserves_stdout_marker():
    assert prepare_output_file("-") == "-"


def test_open_output_file_uses_stdout_for_dash(monkeypatch):
    stdout = io.StringIO()
    monkeypatch.setattr("sys.stdout", stdout)

    with open_output_file("-") as (output_file, output_handle):
        output_handle.write("hello stdout")

    assert output_file == "-"
    assert stdout.getvalue() == "hello stdout"
