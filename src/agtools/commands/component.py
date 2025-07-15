#!/usr/bin/env python3

import re

from agtools.core.graph import UnitigGraph
from agtools.log_config import logger


def _write_component_graph(component_segments, gfa_file, output_path):

    output_file = f"{output_path}/component_graph.gfa"

    with open(gfa_file, "r") as gfa, open(output_file, "w") as filtered_gfa:
        for line in gfa:
            if line.startswith("S"):
                parts = line.strip().split("\t")
                seg_id = parts[1]
                if seg_id in component_segments:
                    filtered_gfa.write(line)
            elif line.startswith("L") or line.startswith("J"):
                parts = line.strip().split("\t")
                from_seg, to_seg = parts[1], parts[3]
                if from_seg in component_segments and to_seg in component_segments:
                    filtered_gfa.write(line)
            elif line.startswith("C"):
                parts = line.strip().split("\t")
                container_seg, contained_seg = parts[1], parts[3]
                if (
                    container_seg in component_segments
                    and contained_seg in component_segments
                ):
                    filtered_gfa.write(line)
            elif line.startswith("P"):
                parts = line.strip().split("\t")
                seg_ids = parts[2].split(",")
                if all(seg_id in component_segments for seg_id in seg_ids):
                    filtered_gfa.write(line)
            elif line.startswith("W"):
                parts = line.strip().split("\t")
                seg_ids = re.split(r"[><]", parts[-1])
                if all(seg_id in component_segments for seg_id in seg_ids):
                    filtered_gfa.write(line)
            else:
                filtered_gfa.write(line)

    return output_file


def component(gfa_file, segment, output):

    ug = UnitigGraph.from_gfa(gfa_file)

    connected_components = ug.graph.components()

    segment_id = ug.segment_names_rev[segment]

    component_segments = []

    for component in connected_components:
        if segment_id in component:
            component_segments = [ug.segment_names[node_id] for node_id in component]
            break

    output_file = _write_component_graph(component_segments, gfa_file, output)

    return output_file
