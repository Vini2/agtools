#!/usr/bin/env python3

import logging
from collections import OrderedDict
from pathlib import Path
from typing import Mapping, Optional

import click

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "0.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Alpha"

# Setup logger
# ---------------------------------------------------

logger = logging.getLogger(f"agtools {__version__}")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
consoleHeader = logging.StreamHandler()
consoleHeader.setFormatter(formatter)
consoleHeader.setLevel(logging.INFO)
logger.addHandler(consoleHeader)

class OrderedGroup(click.Group):
    """custom group class to ensure help function returns commands in desired order.
    class is adapted from Максим Стукало's answer to
    https://stackoverflow.com/questions/47972638/how-can-i-define-the-order-of-click-sub-commands-in-help
    """

    def __init__(
        self,
        name: Optional[str] = None,
        commands: Optional[Mapping[str, click.Command]] = None,
        **kwargs,
    ):
        super().__init__(name, commands, **kwargs)
        #: the registered subcommands by their exported names.
        self.commands = commands or OrderedDict()

    def list_commands(self, ctx: click.Context) -> Mapping[str, click.Command]:
        return self.commands


@click.group(
    cls=OrderedGroup, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(__version__, "-v", "--version", is_flag=True)
def main():
    """agtools: Tools for manipulating assembly graphs"""
    pass

_graph = click.option(
    "--graph",
    help="path to the assembly graph file",
    type=click.Path(exists=True),
    required=True,
)
_output = click.option(
    "--output",
    help="path to the output folder",
    type=click.Path(dir_okay=True, writable=True, readable=True),
    required=True,
)
_prefix = click.option(
    "--prefix",
    help="prefix for the output file",
    type=str,
    default="",
    required=False,
)


_click_command_opts = dict(
    no_args_is_help=True, context_settings={"show_default": True}
)


@main.command(**_click_command_opts)
@_graph
@_output
def gfa2fasta(
    graph,
    output
):
    """Get segments in FASTA format"""
    print("Running gfa2fasta")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--ksize",
    help="k-mer size used for the assembly",
    type=int,
    default=141,
    show_default=True,
    required=True,
)
@_output
def fastg2gfa(
    graph,
    ksize,
    output
):
    """Convert FASTG file to GFA"""
    print("Running fastg2gfa")


@main.command(**_click_command_opts)
@_graph
@_output
def rename(
    graph,
    output
):
    """Rename segments in a GFA file"""
    print("Renaming segments")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--segment",
    help="segment ID",
    type=str,
    show_default=True,
    required=True,
)
@_output
def subgraph(
    graph,
    segment,
    output
):
    """Extract a subgraph containing a given segment"""
    print("Extracting a subgraph given a segment")