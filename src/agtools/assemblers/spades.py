#!/usr/bin/env python3
import re
from collections import defaultdict
from agtools.core.graph import UnitigGraph, ContigGraph

from bidict import bidict
from igraph import Graph


def _get_segment_paths(contig_paths):
    paths = {}
    segment_contigs = {}
    node_count = 0

    id_map = bidict()       # id → contig_num
    contig_names = bidict() # id → contig_name

    current_contig_num = ""

    with open(contig_paths) as file:
        name = file.readline().strip()
        path = file.readline().strip()

        while name != "" and path != "":
            while ";" in path:
                path = path[:-2] + "," + file.readline()

            start = "NODE_"
            end = "_length_"
            contig_num = str(int(re.search("%s(.*)%s" % (start, end), name).group(1)))

            segments = path.rstrip().split(",")

            if current_contig_num != contig_num:
                id_map[node_count] = int(contig_num)
                contig_names[node_count] = name.strip()
                current_contig_num = contig_num
                node_count += 1

            if contig_num not in paths:
                paths[contig_num] = segments

            for segment in segments:
                if segment not in segment_contigs:
                    segment_contigs[segment] = set([contig_num])
                else:
                    segment_contigs[segment].add(contig_num)

            name = file.readline().strip()
            path = file.readline().strip()

    return paths, segment_contigs, node_count, id_map, contig_names


def _get_graph_edges(
    graph_file, contigs_map, contigs_map_rev, paths, segment_contigs
):
    links = []
    links_map = defaultdict(set)

    # Get links from assembly_graph_with_scaffolds.gfa
    with open(graph_file) as file:
        line = file.readline()

        while line != "":
            # Identify lines with link information
            if "L" in line:
                strings = line.split("\t")
                f1, f2 = strings[1] + strings[2], strings[3] + strings[4]
                links_map[f1].add(f2)
                links_map[f2].add(f1)
                links.append(strings[1] + strings[2] + " " + strings[3] + strings[4])
            line = file.readline()

    # Create list of edges
    edge_list = []

    for i in range(len(paths)):
        segments = paths[str(contigs_map[i])]

        new_links = []

        for segment in segments:
            my_segment = segment

            my_segment_rev = ""

            if my_segment.endswith("+"):
                my_segment_rev = my_segment[:-1] + "-"
            else:
                my_segment_rev = my_segment[:-1] + "+"

            if segment in links_map:
                new_links.extend(list(links_map[segment]))

            if my_segment_rev in links_map:
                new_links.extend(list(links_map[my_segment_rev]))

        if my_segment in segment_contigs:
            for contig in segment_contigs[my_segment]:
                if i != contigs_map_rev[int(contig)]:
                    # Add edge to list of edges
                    edge_list.append((i, contigs_map_rev[int(contig)]))

        if my_segment_rev in segment_contigs:
            for contig in segment_contigs[my_segment_rev]:
                if i != contigs_map_rev[int(contig)]:
                    # Add edge to list of edges
                    edge_list.append((i, contigs_map_rev[int(contig)]))

        for new_link in new_links:
            if new_link in segment_contigs:
                for contig in segment_contigs[new_link]:
                    if i != contigs_map_rev[int(contig)]:
                        # Add edge to list of edges
                        edge_list.append((i, contigs_map_rev[int(contig)]))

    return edge_list


def get_contig_graph(graph_file, contig_paths_file) -> ContigGraph:
    # Get paths, segments, links and contigs of the assembly graph
    (
        contig_paths,
        segment_contigs,
        node_count,
        contigs_map,
        contig_names,
    ) = _get_segment_paths(contig_paths_file)

    # Create graph
    graph = Graph()

    # Add vertices
    graph.add_vertices(node_count)

    # Name vertices with contig identifiers
    for i in range(node_count):
        graph.vs[i]["id"] = i
        graph.vs[i]["label"] = contig_names[i]

    # Get list of edges
    edge_list = _get_graph_edges(
        graph_file=graph_file,
        contigs_map=contigs_map,
        contigs_map_rev=contigs_map.inverse,
        paths=contig_paths,
        segment_contigs=segment_contigs,
    )

    # Add edges to the graph
    graph.add_edges(edge_list)

    # Simplify the graph
    graph.simplify(multiple=True, loops=False, combine_edges=None)

    contig_graph = ContigGraph(
        graph = graph,
        path = graph_file,
        contig_ids = contigs_map,
        contig_names = contig_names,
        graph_to_contig_map = None,
    )

    return contig_graph


def get_unitig_graph(graph_file) -> UnitigGraph:
    ug = UnitigGraph.from_gfa(graph_file)
    return ug