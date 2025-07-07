#!/usr/bin/env python3

import sys
from collections import OrderedDict
from pathlib import Path
from typing import Mapping, Optional

import click
from loguru import logger

from agtools import commands

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

# Remove the default logger configuration
logger.remove()

# Console logging (INFO level and above)
logger.add(sink=sys.stdout, level="INFO")

# File logging (DEBUG level and above)
logger.add(sink="agtools.log", level="DEBUG")


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
    "-g",
    help="path to the assembly graph file",
    type=click.Path(exists=True),
    required=True,
)
_output = click.option(
    "--output",
    "-o",
    help="path to the output folder",
    type=click.Path(dir_okay=True, writable=True, readable=True),
    required=True,
)


_click_command_opts = dict(
    no_args_is_help=True, context_settings={"show_default": True}
)


@main.command(**_click_command_opts)
@_graph
@_output
def stats(graph, output):
    """Compute statistics about the graph"""
    print("Running stats")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--prefix",
    "-p",
    help="prefix for the output file",
    type=str,
    default="",
    required=False,
)
@_output
def rename(graph, prefix, output):
    """Rename segments in a GFA file"""

    logger.info(f"Renaming segments in graph file {graph}")
    logger.info(f"Prefix used is {prefix}")

    output_file = commands.rename(graph, prefix, output)

    logger.info(f"Renamed graph file is {output_file}")


@main.command(**_click_command_opts)
@_graph
@_output
def merge(graph, output):
    """Merge two or more GFA files"""
    print("Merging GFA files")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--min-length",
    "-l",
    help="minimum length of segments to filter",
    type=int,
    default=141,
    show_default=True,
    required=True,
)
@_output
def filter(graph, length, output):
    """Filter segments from GFA file"""
    print("Running filter")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--segment",
    "-s",
    help="segment ID",
    type=str,
    show_default=True,
    required=True,
)
@_output
def component(graph, segment, output):
    """Extract a component containing a given segment"""
    print("Extracting a component given a segment")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--ksize",
    "-k",
    help="k-mer size used for the assembly",
    type=int,
    default=141,
    show_default=True,
    required=True,
)
@_output
def fastg2gfa(graph, ksize, output):
    """Convert FASTG file to GFA format"""

    logger.info(f"Converting FASTG file {graph} to GFA format")
    logger.info(f"k-mer size {ksize} will be used as the overlap")

    gfa_path = commands.fastg2gfa(graph, ksize, output)

    logger.info(f"GFA file written to {gfa_path} with fixed overlap: {ksize}M")


@main.command(**_click_command_opts)
@_graph
@_output
def gfa2fastg(graph, output):
    """Convert GFA file to FASTG format"""

    logger.info(f"Converting GFA file {graph} to FASTG format")

    fastg_path, overlap_value = commands.gfa2fastg(graph, output)

    logger.info(f"The detected overlap value is {overlap_value}")
    logger.info(f"FASTG file written to {fastg_path}")


@main.command(**_click_command_opts)
@_graph
@_output
def gfa2dot(graph, ksize, output):
    """Convert GFA file to DOT format (Graphviz)"""
    print("Running gfa2dot")


@main.command(**_click_command_opts)
@_graph
@_output
def gfa2fasta(graph, output):
    """Get segments in FASTA format"""

    logger.info(f"Extracting segment sequences from {graph} file in to FASTA format")

    fasta_path = commands.gfa2fasta(graph, output)

    logger.info(f"FASTA file written to {fasta_path}")


@main.command(**_click_command_opts)
@_graph
@_output
def gfa2adj(graph, output):
    """Get adjacency matrix of the assembly graph"""
    print("Running gfa2adj")
