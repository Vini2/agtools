#!/usr/bin/env python3

import io
from collections.abc import Callable

from igraph import Graph

from agtools.core.contig_graph import ContigGraph
from agtools.core.fasta_parser import FastaParser
from agtools.core.unitig_graph import UnitigGraph


def get_segments(graph_file: str) -> dict[str, int]:
    """
    Parse a GFA file to extract segment names and their corresponding IDs.
    """

    segment_name_to_id = {}
    segment_names = []

    with io.open(graph_file, mode="r", buffering=1024 * 1024) as f:
        while True:
            line = f.readline()
            if not line:
                break

            tag = line[0]
            if tag == "S":  # Segment line
                parts = line.rstrip().split("\t")
                seg_name = parts[1]
                seg_id = len(segment_names)
                segment_name_to_id[seg_name] = seg_id
                segment_names.append(seg_name)

    return segment_name_to_id


def get_graph_edges(
    graph_file: str, segment_contigs: dict, segment_name_to_id: dict
) -> tuple[list[tuple[int, int]], list[int], int]:
    """
    Construct edges between contigs based on shared segment links in the GFA file.
    """

    lcount = 0
    self_loops = set()
    edge_list = set()

    with io.open(graph_file, mode="r", buffering=1024 * 1024) as file:
        line = file.readline()

        while line != "":
            # Identify lines with link information
            if "L" in line:
                lcount += 1
                strings = line.split("\t")
                source = segment_name_to_id[strings[1]]
                target = segment_name_to_id[strings[3]]

                source_contigs = segment_contigs.get(source)
                target_contigs = segment_contigs.get(target)

                if source_contigs and target_contigs:
                    for source_contig in source_contigs:
                        for target_contig in target_contigs:
                            if (
                                source_contig != target_contig
                                and (source_contig, target_contig) not in edge_list
                            ):
                                edge_list.add((source_contig, target_contig))
                            else:
                                self_loops.add(source_contig)

            line = file.readline()

    return list(edge_list), list(self_loops), lcount


def build_contig_graph(
    graph_file: str,
    contigs_file: str,
    contig_paths_file: str,
    segment_path_parser: Callable[[str, dict], tuple[dict, list, dict]],
) -> ContigGraph:
    """
    Build a contig-level graph from a GFA file and parser-specific contig paths file.
    """

    segment_name_to_id = get_segments(graph_file)

    segment_contigs, contig_names, contig_name_to_id = segment_path_parser(
        contig_paths_file, segment_name_to_id
    )

    graph = Graph()
    graph.add_vertices(len(contig_names))
    graph.vs["label"] = contig_names

    edge_list, self_loops, lcount = get_graph_edges(
        graph_file=graph_file,
        segment_contigs=segment_contigs,
        segment_name_to_id=segment_name_to_id,
    )

    graph.add_edges(edge_list)
    graph.simplify(multiple=True, loops=False, combine_edges=None)

    parser = FastaParser(contigs_file)

    return ContigGraph(
        graph=graph,
        vcount=graph.vcount(),
        lcount=lcount,
        ecount=graph.ecount(),
        file_path=graph_file,
        contig_names=contig_names,
        contig_name_to_id=contig_name_to_id,
        contig_parser=parser,
        contig_descriptions=None,
        graph_to_contig_map=None,
        self_loops=self_loops,
    )


def get_unitig_graph(graph_file: str) -> UnitigGraph:
    """
    Build a unitig-level assembly graph from a GFA file.
    """

    return UnitigGraph.from_gfa(graph_file)
