#!/usr/bin/env python3

import re

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = "1.0.1"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


def _remap_element(element_id: str, element_map: dict) -> str:
    """
    Remap an element ID using the provided mapping.

    Parameters
    ----------
    element_id : str
        The original element ID to be remapped.

    element_map : dict
        Dictionary mapping original element IDs to new element IDs.

    Returns
    -------
    str
        The remapped element ID if found in the mapping, otherwise the original ID.
    """

    return element_map.get(element_id, element_id)


def _build_element_maps(input_gfa: str, prefix: str) -> tuple:
    """
    Create a mapping of element IDs from an input GFA file, applying
    a prefix to each element ID. Used for segments, paths and walks.

    Parameters
    ----------
    input_gfa : str
        Path to the input GFA file.

    prefix : str
        Prefix to prepend to each element ID.

    Returns
    -------
    A tuple containing:
        - segment_map : dict[str, str]
            A dictionary mapping original segment IDs to prefixed segment IDs.
        - path_map : dict[str, str]
            A dictionary mapping original path IDs to prefixed path IDs.
        - walk_map : dict[str, str]
            A dictionary mapping original walk IDs to prefixed walk IDs.
    """

    segment_map = {}
    path_map = {}
    walk_map = {}

    # Build map of old_id -> new_id
    with open(input_gfa, "r") as infile:
        for line in infile:
            if line.startswith("S"):
                parts = line.strip().split("\t")
                old_id = parts[1]
                new_id = f"{prefix}_{old_id}"
                segment_map[old_id] = new_id

            elif line.startswith("P"):
                parts = line.strip().split("\t")
                old_id = parts[1]
                new_id = f"{prefix}_{old_id}"
                path_map[old_id] = new_id

            elif line.startswith("W"):
                parts = line.strip().split("\t")
                old_id = parts[1]
                new_id = f"{prefix}_{old_id}"
                walk_map[old_id] = new_id

    return segment_map, path_map, walk_map


def _write_renamed_file(
    input_gfa: str, segment_map: dict, path_map: dict, walk_map: dict, output_path: str
) -> str:
    """
    Write a new GFA file with segment IDs renamed based on the provided segment map.

    Parameters
    ----------
    input_gfa : str
        Path to the original GFA file.

    segment_map : dict
        Mapping of old segment IDs to new IDs.

    path_map : dict
        Mapping of old path IDs to new IDs.

    walk_map : dict
        Mapping of old walk IDs to new IDs.

    output_path : str
        Directory path where the renamed GFA file will be saved.

    Returns
    -------
    str
        Path to the renamed GFA file.
    """

    output_file = f"{output_path}/renamed_graph.gfa"

    # Rewrite file with renamed segment IDs
    with open(input_gfa, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            parts = line.strip().split("\t")
            if not parts:
                continue

            tag = parts[0]

            if tag == "S":
                parts[1] = _remap_element(parts[1], segment_map)
                outfile.write("\t".join(parts) + "\n")

            elif tag == "L" or tag == "J" or tag == "C":
                parts[1] = _remap_element(parts[1], segment_map)
                parts[3] = _remap_element(parts[3], segment_map)
                outfile.write("\t".join(parts) + "\n")

            elif tag == "P":
                parts[1] = _remap_element(parts[1], path_map)
                path_path = parts[2]
                segments = re.split(r"([,;])", path_path)
                segments = [
                    _remap_element(s[:-1], segment_map) + s[-1] for s in segments
                ]
                parts[2] = "".join(segments)
                outfile.write("\t".join(parts) + "\n")

            elif tag == "W":
                parts[1] = _remap_element(parts[1], walk_map)
                walk_path = parts[-1]
                segments = re.split(r"([><])", walk_path)
                segments = [_remap_element(s, segment_map) for s in segments]
                parts[-1] = "".join(segments)
                outfile.write("\t".join(parts) + "\n")

            else:
                outfile.write(line)

    return output_file


def rename(gfa_file: str, prefix: str, output_path: str) -> str:
    """
    Rename segment IDs in a GFA file by applying a prefix and save the modified file.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.

    prefix : str
        Prefix to prepend to each segment ID.

    output_path : str
        Directory path where the renamed GFA file will be saved.

    Returns
    -------
    str
        Path to the renamed GFA file.
    """

    segment_map, path_map, walk_map = _build_element_maps(gfa_file, prefix)
    output_file = _write_renamed_file(
        gfa_file, segment_map, path_map, walk_map, output_path
    )
    return output_file
