#!/usr/bin/env python3

import io

from agtools import __version__
from agtools.assemblers.flye import get_contig_links, get_contig_paths
from agtools.commands._format_checks import validate_gfa_input
from agtools.commands._output import open_output_file
from agtools.core.fasta_parser import FastaParser
from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = __version__
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def _read_assembly_info(info_file: str) -> dict:
    """
    Parse the metadata table that Flye writes alongside its assembly.

    The ``assembly_info.txt`` file is a tab-separated table whose first line is
    a header starting with ``#seq_name``. Only the columns needed to annotate
    the contig-level graph are retained.

    Parameters
    ----------
    info_file : str
        Path to the Flye ``assembly_info.txt`` file.

    Returns
    -------
    dict[str, dict]
        Mapping from contig name to a dictionary with the ``length``,
        ``coverage`` and ``circular`` entries for that contig.
    """

    contig_info = {}

    with io.open(info_file, mode="r", buffering=1024 * 1024) as file:
        for line in file:
            if line.startswith("#") or line.startswith("seq_name"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 4:
                continue

            contig_info[parts[0]] = {
                "length": int(parts[1]),
                "coverage": int(float(parts[2])),
                "circular": parts[3].strip().upper() == "Y",
            }

    return contig_info


def _write_contig_gfa(
    contig_names: list,
    parser: FastaParser,
    contig_info: dict,
    links: list,
    output_path: str,
) -> str:
    """
    Write a contig-level GFA file.

    Each contig becomes one ``S`` record carrying its sequence from the FASTA
    file, and each derived adjacency becomes one ``L`` record. Sequence lengths
    are always tagged with ``LN:i:``. When assembly metadata is available, read
    coverage is tagged with ``dp:i:`` and the topology reported by the
    assembler is tagged with ``TP:Z:circular`` or ``TP:Z:linear``.

    A self-link is a genuine adjacency of the assembly graph and is written
    whenever one is found, but it does not by itself mean the contig is
    circular: a repeat contig whose two ends are adjacent in the graph carries
    one too. The ``TP:Z:`` tag is the authoritative circularity call.

    Parameters
    ----------
    contig_names : list[str]
        Contig names to write, in output order.
    parser : FastaParser
        Parser for the FASTA file holding the contig sequences.
    contig_info : dict
        Contig metadata as returned by :func:`_read_assembly_info`.
    links : list[tuple[str, str, str, str, str]]
        Contig-level links.
    output_path : str
        Path where the contig-level GFA file will be written.

    Returns
    -------
    str
        Full path to the written GFA file.
    """

    with open_output_file(output_path) as (output_file, output_handle):
        output_handle.write("H\tVN:Z:1.0\n")

        for contig_name in contig_names:
            sequence = str(parser.get_sequence(contig_name))

            fields = ["S", contig_name, sequence, f"LN:i:{len(sequence)}"]

            info = contig_info.get(contig_name)
            if info is not None:
                fields.append(f"dp:i:{info['coverage']}")
                fields.append("TP:Z:circular" if info["circular"] else "TP:Z:linear")

            output_handle.write("\t".join(fields) + "\n")

        for from_contig, from_orient, to_contig, to_orient, overlap in links:
            output_handle.write(
                "\t".join(
                    ["L", from_contig, from_orient, to_contig, to_orient, overlap]
                )
                + "\n"
            )

    return output_file


def flye2contig(
    gfa_file: str,
    fasta: str,
    info: str | None,
    output_path: str,
) -> str:
    """
    Convert a Flye edge-level assembly graph into a contig-level GFA file.

    Flye writes an assembly graph whose ``S`` records are graph edges
    (``edge_X``) rather than the contigs (``contig_X``) of ``assembly.fasta``.
    The contigger merges chains of edges into single contigs and extends them
    into flanking repeats, so the graph holds more segments than the assembly
    holds contigs, under different names.

    This function reads the ``P`` records that record which oriented edges each
    contig walks through, takes every contig sequence verbatim from the FASTA
    file, and lifts the edge-level ``L`` records to contig-level links. Contigs
    that appear as a path in the graph but not in the FASTA file were discarded
    by Flye and are dropped from the output.

    .. note::
       Only Flye assemblies are supported for now. The graph must be a raw Flye
       ``assembly_graph.gfa``, which carries one ``P`` record per contig.

    Parameters
    ----------
    gfa_file : str
        Path to the Flye edge-level GFA file (``assembly_graph.gfa``).
    fasta : str
        Path to the FASTA file with the final contig sequences
        (``assembly.fasta``).
    info : str or None
        Optional path to the Flye assembly metadata table
        (``assembly_info.txt``), used to tag coverage and topology and to check
        the contig lengths of the output.
    output_path : str
        Path where the contig-level GFA file will be saved.

    Returns
    -------
    str
        Full path to the contig-level GFA file.

    Raises
    ------
    ValueError
        If the input is not a GFA file, or if the GFA file holds no contig
        paths.
    """

    validate_gfa_input(gfa_file, "flye2contig")

    contig_paths = get_contig_paths(gfa_file)

    if not contig_paths:
        message = (
            f"No contig paths (P records) were found in {gfa_file}. "
            "The flye2contig subcommand expects a raw Flye assembly graph "
            "that records the edges walked by each contig."
        )
        logger.error(message)
        raise ValueError(message)

    logger.info(f"Found {len(contig_paths)} contig paths in the assembly graph")

    parser = FastaParser(fasta)

    contig_names = [name for name in parser.index if name in contig_paths]

    missing_paths = [name for name in parser.index if name not in contig_paths]
    if missing_paths:
        logger.warning(
            f"{len(missing_paths)} contigs in {fasta} have no path in the "
            "assembly graph and are written without any links"
        )
        contig_names.extend(missing_paths)

    dropped = len(contig_paths) - (len(contig_names) - len(missing_paths))
    if dropped > 0:
        logger.warning(
            f"Dropping {dropped} contig paths that are not present in {fasta}. "
            "These were discarded by Flye when writing the final contigs."
        )

    kept_paths = {
        name: contig_paths[name] for name in contig_names if name in contig_paths
    }

    contig_info = _read_assembly_info(info) if info else {}

    if contig_info:
        for contig_name in contig_names:
            expected = contig_info.get(contig_name, {}).get("length")
            observed = len(parser.get_sequence(contig_name))
            if expected is not None and expected != observed:
                logger.warning(
                    f"Length of {contig_name} is {observed} bp in {fasta} but "
                    f"{expected} bp in {info}"
                )

    links = get_contig_links(gfa_file, kept_paths)

    if contig_info:
        self_linked = {
            from_contig
            for from_contig, from_orient, to_contig, to_orient, _ in links
            if from_contig == to_contig and from_orient == to_orient
        }
        circular = {
            name for name in contig_names if contig_info.get(name, {}).get("circular")
        }
        logger.info(f"{len(circular)} contigs are reported as circular by Flye")

        if self_linked - circular:
            logger.warning(
                f"{len(self_linked - circular)} contigs have a self-link in the "
                f"assembly graph but are not circular according to {info}. "
                "Use the TP:Z: tag rather than the presence of a self-link to "
                "decide whether a contig is circular."
            )

    logger.info(f"Writing {len(contig_names)} contigs and {len(links)} links")

    output_file = _write_contig_gfa(
        contig_names, parser, contig_info, links, output_path
    )

    return output_file
