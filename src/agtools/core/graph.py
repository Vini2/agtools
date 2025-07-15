#!/usr/bin/env python3

import re
from collections import defaultdict

from bidict import bidict
from Bio.Seq import Seq
from igraph import Graph


class AssemblyGraph:
    def __init__(self):
        self.graph = Graph(directed=False)
        self.path = None
        self.oriented_links = defaultdict(lambda: defaultdict(list))
        self.link_overlap = {}
        self.segment_names = None  # node_id → segment_id
        self.segment_names_rev = None  # segment_id → node_id
        self.segment_sequences = {}  # segment_id → sequence
        self.segment_lengths = {}  # segment_id → length
        self.self_loops = []

    @classmethod
    def from_gfa(cls, path: str) -> "AssemblyGraph":
        ag = cls()
        node_count = 0
        links = []
        segment_names = {}

        ag.path = path

        with open(path) as f:
            for line in f:

                if line.startswith("S"):
                    parts = line.strip().split("\t")
                    seg_id = parts[1]
                    seq = parts[2]
                    segment_names[node_count] = seg_id
                    ag.segment_sequences[seg_id] = Seq(seq)
                    ag.segment_lengths[seg_id] = len(seq)
                    node_count += 1
                elif line.startswith("L"):
                    parts = line.strip().split("\t")
                    from_seg, from_orient = parts[1], parts[2]
                    to_seg, to_orient = parts[3], parts[4]
                    overlap = int(parts[5][:-1])  # Remove trailing M

                    links.append((from_seg, to_seg))
                    ag._add_oriented_links(
                        from_seg, to_seg, from_orient, to_orient, overlap
                    )

        ag.segment_names = bidict(segment_names)
        ag.segment_names_rev = ag.segment_names.inverse
        ag.graph.add_vertices(node_count)

        for i in range(node_count):
            seg_id = ag.segment_names[i]
            ag.graph.vs[i]["id"] = i
            ag.graph.vs[i]["name"] = seg_id
            ag.graph.vs[i]["label"] = f"{seg_id}\nID:{i}"

        edge_list, ag.self_loops = ag._get_graph_edges(links)
        ag.graph.add_edges(edge_list)
        ag.graph.simplify(multiple=True, loops=False)
        return ag

    def _add_oriented_links(self, from_seg, to_seg, from_orient, to_orient, overlap):
        key1 = f"{from_seg}{from_orient}"
        key2 = f"{to_seg}{to_orient}"
        self.oriented_links[from_seg][to_seg].append((from_orient, to_orient))
        self.link_overlap[(key1, key2)] = overlap

        # Add symmetric reverse
        rev1 = "+" if from_orient == "-" else "-"
        rev2 = "+" if to_orient == "-" else "-"
        self.oriented_links[to_seg][from_seg].append((rev2, rev1))
        self.link_overlap[(f"{to_seg}{rev2}", f"{from_seg}{rev1}")] = overlap

        print(self.oriented_links)

    def _get_graph_edges(self, links):
        edges = []
        loops = []
        for from_edge, to_edge in links:
            if from_edge == to_edge:
                loops.append(from_edge)
            else:
                src = self.segment_names_rev[from_edge]
                tgt = self.segment_names_rev[to_edge]
                edges.append((src, tgt))
        return edges, loops

    def get_neighbors(self, seg_id: str) -> list:
        vid = self.segment_names_rev[seg_id]
        neighbor_ids = self.graph.neighbors(vid)
        return [self.segment_names[nid] for nid in neighbor_ids]

    def filter_segments(self, min_length: int) -> "AssemblyGraph":
        keep_segs = {
            s for s, seq in self.segment_sequences.items() if len(seq) >= min_length
        }
        keep_ids = [self.segment_names_rev[s] for s in keep_segs]
        subgraph = self.graph.subgraph(keep_ids)

        new_ag = AssemblyGraph()
        new_ag.graph = subgraph
        new_ag.segment_names = bidict(
            {i: self.segment_names[v.index] for i, v in enumerate(subgraph.vs)}
        )
        new_ag.segment_names_rev = new_ag.segment_names.inverse
        new_ag.segment_sequences = {s: self.segment_sequences[s] for s in keep_segs}
        new_ag.segment_lengths = {
            s: len(seq) for s, seq in new_ag.segment_sequences.items()
        }
        return new_ag


def parse_fastg(fastg_file):
    segments = {}
    edges = {}

    with open(fastg_file, "r") as f:
        current_node = None
        sequence = []

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_node and sequence:
                    segments[current_node] = "".join(sequence)
                    sequence = []

                header = line[1:]
                parts = header.split(":")
                node = parts[0].strip("'")
                neighbors = []
                if len(parts) > 1:
                    neighbors = re.split(r"[,\s]+", parts[1])
                current_node = node
                edges[node] = neighbors
            else:
                sequence.append(line)

        if current_node and sequence:
            segments[current_node] = "".join(sequence)

    return segments, edges
