#!/usr/bin/env python3

import re

from agtools import __version__
from agtools.commands._output import open_output_file

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]
__license__ = "MIT"
__version__ = __version__
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


DEFAULT_SEPARATOR = "_"


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


def _check_collisions(element_map: dict, element_type: str) -> None:
    """
    Check that renaming does not produce IDs that clash with existing ones.

    A collision occurs when a newly generated ID is already used by a
    different element of the same type in the input file, which would make
    the two elements indistinguishable in the renamed file.

    Parameters
    ----------
    element_map : dict
        Mapping of original element IDs to new element IDs.

    element_type : str
        Name of the element type being checked (e.g. "segment"), used in the
        error message.

    Raises
    ------
    ValueError
        If any new ID clashes with an existing ID, or if two different
        original IDs map onto the same new ID.
    """

    existing_ids = set(element_map)

    clashes = sorted(
        old_id
        for old_id, new_id in element_map.items()
        if new_id != old_id and new_id in existing_ids
    )

    if clashes:
        raise ValueError(
            f"Renaming would produce duplicate {element_type} IDs. "
            f"The following {element_type} IDs already exist in the input file "
            f"after renaming: {', '.join(element_map[old_id] for old_id in clashes)}. "
            "Please choose a different prefix or separator."
        )

    seen = {}

    for old_id, new_id in element_map.items():
        if new_id in seen:
            raise ValueError(
                f"Renaming would produce duplicate {element_type} IDs. "
                f"Both '{seen[new_id]}' and '{old_id}' map to '{new_id}'. "
                "Please choose a different prefix or separator."
            )

        seen[new_id] = old_id


def _build_element_maps(
    input_gfa: str, prefix: str, separator: str = DEFAULT_SEPARATOR
) -> tuple:
    """
    Create a mapping of element IDs from an input GFA file, applying
    a prefix to each element ID. Used for segments, paths and walks.

    An empty prefix leaves all IDs unchanged.

    Parameters
    ----------
    input_gfa : str
        Path to the input GFA file.

    prefix : str
        Prefix to prepend to each element ID.

    separator : str
        String placed between the prefix and the original ID.

    Returns
    -------
    A tuple containing:
        - segment_map : dict[str, str]
            A dictionary mapping original segment IDs to prefixed segment IDs.
        - path_map : dict[str, str]
            A dictionary mapping original path IDs to prefixed path IDs.
        - walk_map : dict[str, str]
            A dictionary mapping original walk IDs to prefixed walk IDs.

    Raises
    ------
    ValueError
        If renaming would produce duplicate segment, path or walk IDs.
    """

    segment_map = {}
    path_map = {}
    walk_map = {}

    maps_by_tag = {"S": segment_map, "P": path_map, "W": walk_map}

    # An empty prefix is a no-op, so no separator is prepended
    new_prefix = f"{prefix}{separator}" if prefix else ""

    # Build map of old_id -> new_id
    with open(input_gfa, "r") as infile:
        for line in infile:
            tag = line[:1]

            if tag in maps_by_tag:
                parts = line.strip().split("\t")
                old_id = parts[1]
                maps_by_tag[tag][old_id] = f"{new_prefix}{old_id}"

    _check_collisions(segment_map, "segment")
    _check_collisions(path_map, "path")
    _check_collisions(walk_map, "walk")

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
        Path where the renamed GFA file will be saved.

    Returns
    -------
    str
        Path to the renamed GFA file.
    """

    # Rewrite file with renamed segment IDs
    with (
        open(input_gfa, "r") as infile,
        open_output_file(output_path) as (
            output_file,
            outfile,
        ),
    ):
        for line in infile:
            parts = line.strip().split("\t")
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
                path_segments = re.split(r"([,;])", path_path)
                remapped_path_segments = [
                    _remap_element(s[:-1], segment_map) + s[-1] for s in path_segments
                ]
                parts[2] = "".join(remapped_path_segments)
                outfile.write("\t".join(parts) + "\n")

            elif tag == "W":
                parts[1] = _remap_element(parts[1], walk_map)
                walk_path = parts[-1]
                walk_segments = re.split(r"([><])", walk_path)
                remapped_walk_segments = [
                    _remap_element(s, segment_map) for s in walk_segments
                ]
                parts[-1] = "".join(remapped_walk_segments)
                outfile.write("\t".join(parts) + "\n")

            else:
                outfile.write(line)

    return output_file


def rename(
    gfa_file: str,
    prefix: str,
    output_path: str,
    separator: str = DEFAULT_SEPARATOR,
) -> str:
    """
    Rename segment IDs in a GFA file by applying a prefix and save the modified file.

    Parameters
    ----------
    gfa_file : str
        Path to the input GFA file.

    prefix : str
        Prefix to prepend to each segment ID. An empty prefix leaves all IDs
        unchanged.

    output_path : str
        Path where the renamed GFA file will be saved.

    separator : str
        String placed between the prefix and the original ID. Defaults to "_".

    Returns
    -------
    str
        Path to the renamed GFA file.

    Raises
    ------
    ValueError
        If renaming would produce duplicate segment, path or walk IDs.
    """

    segment_map, path_map, walk_map = _build_element_maps(gfa_file, prefix, separator)
    output_file = _write_renamed_file(
        gfa_file, segment_map, path_map, walk_map, output_path
    )
    return output_file
