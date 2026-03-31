#!/usr/bin/env python3

from agtools.commands._output import prepare_output_file


def test_prepare_output_file_creates_parent_directories_and_expands_user(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("HOME", str(tmp_path))

    output_file = prepare_output_file("~/nested/outputs/result.txt")

    assert output_file == str(tmp_path / "nested" / "outputs" / "result.txt")
    assert (tmp_path / "nested" / "outputs").is_dir()
