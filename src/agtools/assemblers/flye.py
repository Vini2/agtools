#!/usr/bin/env python3

import io
from collections import defaultdict

from agtools.assemblers._contig_graph_base import (
    build_contig_graph as _build_contig_graph,
)
from agtools.assemblers._contig_graph_base import get_unitig_graph as _get_unitig_graph
from agtools.core.contig_graph import ContigGraph
from agtools.core.unitig_graph import UnitigGraph


def _get_segment_paths_and_contig_mapping(
    contig_paths: str, segment_name_to_id: dict
) -> tuple:
    """
    Parse a contig paths file and extract segment-contig relationships.

    Parameters
    ----------
    contig_paths : str
        Path to the contig paths file (e.g. contigs.paths of scaffolds.paths).
    segment_name_to_id : dict[str, int]
        Mapping from segment name to its internal ID.

    Returns
    -------
    tuple
    segment_contigs : dict[str, set[str]]
        Mapping from segment ID to the set of contig numbers it appears in.
    contig_names : list
        List of contig names.
    contig_name_to_id : dict[str, int]
        Mapping from contig name to its internal ID.
    """

    contig_names = []
    contig_name_to_id = dict()

    segment_contigs = defaultdict(set)

    with io.open(contig_paths, mode="r", buffering=1024 * 1024) as file:
        for line in file.readlines():
            if not (line.startswith("#") or line.startswith("seq_name")):
                strings = line.strip().split()

                contig_name = strings[0]
                contig_id = len(contig_names)
                contig_name_to_id[contig_name] = contig_id
                contig_names.append(contig_name)

                path = strings[-1]
                path = path.replace("*", "")

                if path.startswith(","):
                    path = path[1:]

                if path.endswith(","):
                    path = path[:-1]

                segments = path.rstrip().split(",")

                for segment in segments:
                    if segment[0] == "-":
                        segment_contigs[segment_name_to_id[f"edge_{segment[1:]}"]].add(
                            contig_id
                        )
                    else:
                        segment_contigs[segment_name_to_id[f"edge_{segment}"]].add(
                            contig_id
                        )

    return segment_contigs, contig_names, contig_name_to_id


def get_contig_graph(
    graph_file: str, contigs_file: str, contig_paths_file: str
) -> ContigGraph:
    """
    Build a contig-level graph from an assembly GFA file and contig path mappings.

    This function parses contig metadata, links, and path structure to construct an
    undirected graph where each node represents a contig and edges represent linkages
    inferred from shared segments or GFA link data.

    Parameters
    ----------
    graph_file : str
        Path to the GFA file.
    contigs_file : str
        Path to the FASTA file with contig sequences.
    contig_paths_file : str
        Path to the file with segment paths used to build contigs.

    Returns
    -------
    ContigGraph
        An object representing the contig-level graph with node metadata.
    """

    return _build_contig_graph(
        graph_file=graph_file,
        contigs_file=contigs_file,
        contig_paths_file=contig_paths_file,
        segment_path_parser=_get_segment_paths_and_contig_mapping,
    )


def get_unitig_graph(graph_file: str) -> UnitigGraph:
    """
    Build a unitig-level assembly graph from a GFA file.

    Parameters
    ----------
    graph_file : str
        Path to the GFA file.

    Returns
    -------
    UnitigGraph
        Parsed unitig graph object.
    """

    return _get_unitig_graph(graph_file)


def _flip(orientation: str) -> str:
    """Return the opposite GFA orientation sign."""

    return "-" if orientation == "+" else "+"


def _path_tips(steps: list) -> tuple:
    """
    Return the two outward-facing segment ends of a contig path.

    A contig path is a list of oriented segments. Each segment has a 5' end
    (``"L"``) and a 3' end (``"R"``) in its forward orientation. Traversing the
    contig forwards leaves the graph through the outward end of the last
    oriented segment; traversing it backwards leaves through the outward end of
    the first oriented segment.

    Parameters
    ----------
    steps : list[tuple[str, str]]
        Oriented segments of the path as ``(segment_name, orientation)``.

    Returns
    -------
    tuple[tuple[str, str], tuple[str, str]]
        The ``(segment_name, end)`` tips reached when leaving the contig in the
        forward and in the reverse orientation, respectively.
    """

    first_segment, first_orientation = steps[0]
    last_segment, last_orientation = steps[-1]

    forward_tip = (last_segment, "R" if last_orientation == "+" else "L")
    reverse_tip = (first_segment, "L" if first_orientation == "+" else "R")

    return forward_tip, reverse_tip


def _canonical_link(
    from_contig: str,
    from_orientation: str,
    to_contig: str,
    to_orientation: str,
    overlap: str,
) -> tuple:
    """
    Return a canonical representation of a contig-level link.

    A GFA link ``L a oa b ob`` describes the same adjacency as
    ``L b flip(ob) a flip(oa)``. Picking the lexicographically smaller of the
    two makes the pair deduplicable.
    """

    forward = (from_contig, from_orientation, to_contig, to_orientation)
    reverse = (to_contig, _flip(to_orientation), from_contig, _flip(from_orientation))

    return (*min(forward, reverse), overlap)


def get_contig_paths(graph_file: str) -> dict:
    """
    Parse the ``P`` records of a Flye GFA file into contig paths.

    Flye writes one ``P`` record per contig of ``assembly.fasta``, listing the
    oriented graph edges (``S`` records) that the contig walks through. For
    example::

        P	contig_2	edge_11-,edge_10-,edge_9-,edge_2+	*

    Parameters
    ----------
    graph_file : str
        Path to the Flye ``assembly_graph.gfa`` file.

    Returns
    -------
    dict[str, list[tuple[str, str]]]
        Mapping from contig name to its path, given as a list of
        ``(segment_name, orientation)`` pairs.

    Examples
    --------
    >>> get_contig_paths("assembly_graph.gfa")["contig_3"]
    [('edge_3', '+'), ('edge_12', '-')]
    """

    contig_paths = {}

    with io.open(graph_file, mode="r", buffering=1024 * 1024) as file:
        for line in file:
            if not line.startswith("P"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 3:
                continue

            steps = []

            for step in parts[2].split(","):
                step = step.strip()

                if not step or step == "*":
                    continue

                if step[-1] in "+-":
                    steps.append((step[:-1], step[-1]))
                else:
                    steps.append((step, "+"))

            if steps:
                contig_paths[parts[1]] = steps

    return contig_paths


def get_contig_links(graph_file: str, contig_paths: dict) -> list:
    """
    Derive contig-level links from the edge-level links of a Flye GFA file.

    Two contigs are adjacent when an ``L`` record of the edge graph joins the
    outward-facing tip of one contig path to the outward-facing tip of another.
    Edge ends that fall in the interior of a contig path are ignored, because
    those adjacencies are already resolved inside the contig. A contig whose
    two tips are joined to each other yields a self-link, which is how a
    circular contig is represented in GFA.

    Parameters
    ----------
    graph_file : str
        Path to the Flye ``assembly_graph.gfa`` file.
    contig_paths : dict[str, list[tuple[str, str]]]
        Contig paths as returned by :func:`get_contig_paths`, restricted to the
        contigs that should appear in the output.

    Returns
    -------
    list[tuple[str, str, str, str, str]]
        Sorted, deduplicated links as
        ``(from_contig, from_orientation, to_contig, to_orientation, overlap)``.

    Examples
    --------
    >>> get_contig_links("assembly_graph.gfa", {"contig_1": [("edge_1", "+")]})
    [('contig_1', '+', 'contig_1', '+', '0M')]
    """

    tips = defaultdict(list)

    for contig_name, steps in contig_paths.items():
        forward_tip, reverse_tip = _path_tips(steps)
        tips[forward_tip].append((contig_name, "+"))
        tips[reverse_tip].append((contig_name, "-"))

    links = set()

    with io.open(graph_file, mode="r", buffering=1024 * 1024) as file:
        for line in file:
            if not line.startswith("L"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 6:
                continue

            from_segment, from_orientation = parts[1], parts[2]
            to_segment, to_orientation = parts[3], parts[4]
            overlap = parts[5]

            from_end = (from_segment, "R" if from_orientation == "+" else "L")
            to_end = (to_segment, "L" if to_orientation == "+" else "R")

            # An L record is symmetric, so both ends are considered as exits.
            for exit_end, entry_end in ((from_end, to_end), (to_end, from_end)):
                for from_contig, exit_orientation in tips.get(exit_end, ()):
                    for to_contig, entry_orientation in tips.get(entry_end, ()):
                        links.add(
                            _canonical_link(
                                from_contig,
                                exit_orientation,
                                to_contig,
                                _flip(entry_orientation),
                                overlap,
                            )
                        )

    return sorted(links)
