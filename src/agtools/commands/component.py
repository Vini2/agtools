#!/usr/bin/env python3

import re

from agtools.core.unitig_graph import UnitigGraph
from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def _write_component_graph(
    component_segments: set, gfa_file: str, output_path: str
) -> str:
    """
    Write a subgraph of the assembly graph containing only the specified segments.

    This function filters a GFA file and writes a new file containing only:
    - Segments (`S`) that are part of the specified component
    - Links (`L`), jumps (`J`), and containments (`C`) where both endpoints are in the component
    - Paths (`P`) and walks (`W`) that reference only segments in the component
    - Any other lines (headers, comments, etc.) are preserved

    Parameters
    ----------
    component_segments : set of str
        The segment IDs that make up the target component.
    gfa_file : str
        Path to the input GFA file.
    output_path : str
        Directory to write the filtered component GFA file.

    Returns
    -------
    str
        Path to the newly written GFA file containing only the specified component.

    References
    ----------
    The GFA Format Specification
    [https://gfa-spec.github.io/GFA-spec/GFA1.html](https://gfa-spec.github.io/GFA-spec/GFA1.html)
    """

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
                seg_ids = [part[:-1] for part in re.split(r"[,;]", parts[2])]
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


def component(gfa_file: str, segment: str, output_path: str) -> str:
    """
    Extract and write the connected component containing a given segment.

    This function identifies the connected component of the assembly graph that contains
    the given segment. It then writes a filtered GFA file that includes only the segments
    and edges belonging to that component.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.
    segment : str
        Segment ID for which to extract the connected component.
    output : str
        Directory where the filtered GFA file will be saved.

    Returns
    -------
    str
        Path to the component-specific GFA output file.
    """

    ug = UnitigGraph.from_gfa(gfa_file)

    connected_components = ug.graph.components()

    segment_id = ug.segment_name_to_id[segment]

    component_segments = set()

    for component in connected_components:
        if segment_id in component:
            component_segments = set(
                [ug.segment_names[node_id] for node_id in component]
            )
            break

    output_file = _write_component_graph(component_segments, gfa_file, output_path)

    return output_file
