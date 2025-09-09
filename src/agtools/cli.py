#!/usr/bin/env python3

from collections import OrderedDict
from typing import Mapping, Optional

import click

from agtools import commands
from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


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
    """agtools: A Software Framework to Manipulate Assembly Graphs"""
    pass


_graph = click.option(
    "--graph",
    "-g",
    help="path(s) to the assembly graph file(s)",
    type=click.Path(exists=True),
    multiple=True,
    required=True,
)
_output = click.option(
    "--output",
    "-o",
    help="path to the output folder",
    type=click.Path(exists=True, dir_okay=True, writable=True, readable=True),
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

    logger.info(f"Computing statistics of the graph file {graph[0]}")

    output_file = commands.stats(graph[0], output)

    logger.info(f"Computed statistics can be found in {output_file}")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--prefix",
    "-p",
    help="prefix for the graph elements",
    type=str,
    default="",
    required=False,
)
@_output
def rename(graph, prefix, output):
    """Rename segments, paths and walks in a GFA file"""

    logger.info(f"Renaming elements in graph file {graph[0]}")
    logger.info(f"Prefix used is {prefix}")

    output_file = commands.rename(graph[0], prefix, output)

    logger.info(f"Renamed graph file is {output_file}")


@main.command(**_click_command_opts)
@_graph
@_output
def concat(graph, output):
    """Concatenate two or more GFA files"""

    logger.info(f"Concatenating graph files [{', '.join(graph)}]")

    output_file = commands.concat(graph, output)

    logger.info(f"Concatenated graph file is {output_file}")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--min-length",
    "-l",
    help="minimum length of segments to keep",
    type=int,
    default=100,
    show_default=True,
    required=True,
)
@_output
def filter(graph, min_length, output):
    """Filter segments from GFA file"""

    logger.info(f"Filtering segments in graph file {graph[0]}")
    logger.info(f"Minimum length of segments to keep is {min_length} bp")

    filtered_gfa = commands.filter(graph[0], min_length, output)

    logger.info(f"Filtered graph file is {filtered_gfa}")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--fasta",
    "-f",
    help="path to the FASTA file",
    type=click.Path(exists=True),
    required=False,
)
@click.option(
    "--assembler",
    "-a",
    help="assembler name (if assembler used is myloasm)",
    type=str,
    show_default=True,
    required=False,
)
@_output
def clean(graph, fasta, assembler, output):
    """Clean a GFA file based on segments in a FASTA file"""

    logger.info(f"Cleaning the graph file {graph[0]}")
    logger.info(f"Using the FASTA file {fasta}")

    cleaned_gfa = commands.clean(graph[0], fasta, assembler, output)

    logger.info(f"Cleaned graph file is {cleaned_gfa}")


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

    logger.info(
        f"Extracting from file {graph[0]} the component containing segment {segment}"
    )

    gfa_path = commands.component(graph[0], segment, output)

    logger.info(f"GFA file of the component is written to {gfa_path}")


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

    logger.info(f"Converting FASTG file {graph[0]} to GFA format")
    logger.info(f"k-mer size {ksize} will be used as the overlap")

    gfa_path = commands.fastg2gfa(graph[0], ksize, output)

    logger.info(f"GFA file is written to {gfa_path} with fixed overlap: {ksize}M")


@main.command(**_click_command_opts)
@_graph
@_output
def asqg2gfa(graph, output):
    """Convert ASQG file to GFA format"""

    logger.info(f"Converting ASQG file {graph[0]} to GFA format")

    gfa_path = commands.asqg2gfa(graph[0], output)

    logger.info(f"GFA file is written to {gfa_path}")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--abyss",
    "-ab",
    help="use the ABySS DOT format for the output",
    is_flag=True,
    default=False,
    show_default=True,
    required=False,
)
@_output
def gfa2dot(graph, abyss, output):
    """Convert GFA file to DOT format (GraphViz)"""

    logger.info(f"Converting GFA file {graph[0]} to DOT format")

    dot_path = commands.gfa2dot(graph[0], abyss, output)

    logger.info(f"DOT file written to {dot_path}")


@main.command(**_click_command_opts)
@_graph
@_output
def gfa2fasta(graph, output):
    """Get segments in FASTA format"""

    logger.info(f"Extracting segment sequences from {graph[0]} file in to FASTA format")

    fasta_path = commands.gfa2fasta(graph[0], output)

    logger.info(f"FASTA file written to {fasta_path}")


@main.command(**_click_command_opts)
@_graph
@click.option(
    "--delimiter",
    help="delimiter for adjacency file. Supports a comma and a tab.",
    type=click.Choice(["comma", "tab"], case_sensitive=False),
    default="comma",
    show_default=True,
    required=False,
)
@_output
def gfa2adj(graph, delimiter, output):
    """Get adjacency matrix of the assembly graph"""

    logger.info(f"Obtaining the adjacency matrix from {graph[0]}")

    adj_path = commands.gfa2adj(graph[0], delimiter, output)

    logger.info(f"Adjacency matrix is written to {adj_path}")
