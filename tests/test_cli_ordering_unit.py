#!/usr/bin/env python3

from collections import OrderedDict

import click
from click.testing import CliRunner

from agtools.cli import OrderedGroup, main


def test_ordered_group_list_commands_preserves_insertion_order():
    cmd_one = click.Command("one")
    cmd_two = click.Command("two")
    group = OrderedGroup(commands=OrderedDict([("two", cmd_two), ("one", cmd_one)]))

    assert list(group.list_commands(None).keys()) == ["two", "one"]


def test_main_help_lists_commands_in_registered_order():
    runner = CliRunner()
    result = runner.invoke(main, ["--help"], catch_exceptions=False)

    assert result.exit_code == 0
    output = result.output
    assert output.find("\n  stats") < output.find("\n  rename")
    assert output.find("\n  rename") < output.find("\n  concat")
    assert output.find("\n  concat") < output.find("\n  filter")
    assert output.find("\n  fastg2gfa") < output.find("\n  gfa2fastg")
    assert output.find("\n  gfa2fastg") < output.find("\n  asqg2gfa")
    assert output.find("\n  asqg2gfa") < output.find("\n  gfa2asqg")


def test_main_short_help_flag_is_supported():
    runner = CliRunner()
    result = runner.invoke(main, ["-h"], catch_exceptions=False)

    assert result.exit_code == 0
    assert "Usage:" in result.output
