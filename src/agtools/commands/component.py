#!/usr/bin/env python3

from agtools import __version__
from agtools.commands._output import prepare_output_file
from agtools.core.gfa_filter import write_filtered_gfa
from agtools.core.unitig_graph import UnitigGraph
from agtools.log_config import logger

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"

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
        Path to write the filtered component GFA file.

    Returns
    -------
    str
        Path to the newly written GFA file containing only the specified component.

    References
    ----------
    The GFA Format Specification
    [https://gfa-spec.github.io/GFA-spec/GFA1.html](https://gfa-spec.github.io/GFA-spec/GFA1.html)
    """

    output_file = prepare_output_file(output_path)

    def keep_segment(seg_id: str) -> bool:
        return seg_id in component_segments

    write_filtered_gfa(gfa_file, output_file, keep_segment)

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
        Path where the filtered GFA file will be saved.

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
            component_segments = {ug.segment_names[node_id] for node_id in component}
            break

    output_file = _write_component_graph(component_segments, gfa_file, output_path)

    return output_file
