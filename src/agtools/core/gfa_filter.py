#!/usr/bin/env python3

import re
from collections.abc import Callable


def parse_path_segment_ids(path_field: str) -> list[str]:
    """
    Extract segment IDs from a GFA path field.

    Path entries are separated with comma/semicolon and may contain
    orientation suffixes (+/-), which are stripped.
    """

    return [
        segment.rstrip("+-") for segment in re.split(r"[,;]", path_field) if segment
    ]


def parse_walk_segment_ids(walk_field: str) -> list[str]:
    """
    Extract segment IDs from a GFA walk field.

    Walk entries are separated by orientation markers (>/<).
    """

    return [segment for segment in re.split(r"[><]", walk_field) if segment]


def _ensure_newline(line: str) -> str:
    return line if line.endswith("\n") else f"{line}\n"


def _write_parts(output_handle, parts: list[str]) -> None:
    output_handle.write("\t".join(parts) + "\n")


def write_filtered_gfa(
    gfa_file: str,
    output_file: str,
    keep_segment: Callable[[str], bool],
    transform_segment: Callable[[list[str]], list[str]] | None = None,
) -> None:
    """
    Write a filtered GFA file based on a segment-keep predicate.

    Supported tags:
    - S: keep/remove by segment ID and optionally transform the record
    - L/J/C: keep only if both segment IDs are kept
    - P/W: keep only if all referenced segment IDs are kept
    - others: copied unchanged
    """

    with open(gfa_file, "r") as gfa, open(output_file, "w") as filtered_gfa:
        for line in gfa:
            parts = line.rstrip("\n").split("\t")
            if not parts or parts == [""]:
                filtered_gfa.write(_ensure_newline(line))
                continue

            tag = parts[0]

            if tag == "S" and len(parts) > 1:
                if keep_segment(parts[1]):
                    if transform_segment is None:
                        filtered_gfa.write(_ensure_newline(line))
                    else:
                        _write_parts(filtered_gfa, transform_segment(parts))

            elif tag in {"L", "J", "C"} and len(parts) > 3:
                from_seg, to_seg = parts[1], parts[3]
                if keep_segment(from_seg) and keep_segment(to_seg):
                    filtered_gfa.write(_ensure_newline(line))

            elif tag == "P" and len(parts) > 2:
                seg_ids = parse_path_segment_ids(parts[2])
                if all(keep_segment(seg_id) for seg_id in seg_ids):
                    filtered_gfa.write(_ensure_newline(line))

            elif tag == "W" and parts:
                seg_ids = parse_walk_segment_ids(parts[-1])
                if all(keep_segment(seg_id) for seg_id in seg_ids):
                    filtered_gfa.write(_ensure_newline(line))

            else:
                filtered_gfa.write(_ensure_newline(line))
