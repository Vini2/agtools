#!/usr/bin/env python3

from pathlib import Path

from agtools.log_config import logger

_GFA_TAGS = {"H", "S", "L", "J", "C", "P", "W"}
_ASQG_TAGS = {"HT", "VT", "ED"}


def _scan_graph_file(file_path: str) -> tuple[str | None, dict[str, int], str]:
    """Inspect the first non-empty record and count format-specific records."""

    first_format = None
    counts = {"gfa_segments": 0, "fastg_headers": 0, "asqg_segments": 0}

    with open(file_path, "r") as file:
        for line in file:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                continue

            if stripped.startswith(">"):
                if first_format is None:
                    first_format = "FASTG"
                counts["fastg_headers"] += 1
                continue

            tag = stripped.split("\t", 1)[0]

            if first_format is None:
                if tag in _GFA_TAGS:
                    first_format = "GFA"
                elif tag in _ASQG_TAGS:
                    first_format = "ASQG"
                else:
                    first_format = tag

            if tag == "S":
                counts["gfa_segments"] += 1
            elif tag == "VT":
                counts["asqg_segments"] += 1

    return first_format, counts, Path(file_path).name


def _raise_logged_value_error(message: str) -> None:
    logger.error(message)
    raise ValueError(message)


def validate_gfa_input(file_path: str, command_name: str) -> None:
    """Ensure the input file looks like GFA and contains segment records."""

    first_format, counts, file_name = _scan_graph_file(file_path)

    if counts["gfa_segments"] > 0:
        return

    if first_format == "FASTG":
        _raise_logged_value_error(
            f"{file_name} looks like a FASTG file. "
            f"The {command_name} subcommand expects GFA input."
        )

    if first_format == "ASQG":
        _raise_logged_value_error(
            f"{file_name} looks like an ASQG file. "
            f"The {command_name} subcommand expects GFA input."
        )

    _raise_logged_value_error(
        f"No GFA segments were found in {file_name}. "
        f"The {command_name} subcommand expects a GFA file containing S records."
    )


def validate_fastg_input(file_path: str, command_name: str) -> None:
    """Ensure the input file looks like FASTG and contains FASTG headers."""

    first_format, counts, file_name = _scan_graph_file(file_path)

    if counts["fastg_headers"] > 0:
        return

    if first_format == "GFA":
        _raise_logged_value_error(
            f"{file_name} looks like a GFA file. "
            f"The {command_name} subcommand expects FASTG input."
        )

    if first_format == "ASQG":
        _raise_logged_value_error(
            f"{file_name} looks like an ASQG file. "
            f"The {command_name} subcommand expects FASTG input."
        )

    _raise_logged_value_error(
        f"No FASTG segment headers were found in {file_name}. "
        f"The {command_name} subcommand expects a FASTG file containing '>' records."
    )


def validate_asqg_input(file_path: str, command_name: str) -> None:
    """Ensure the input file looks like ASQG and contains VT records."""

    first_format, counts, file_name = _scan_graph_file(file_path)

    if counts["asqg_segments"] > 0:
        return

    if first_format == "GFA":
        _raise_logged_value_error(
            f"{file_name} looks like a GFA file. "
            f"The {command_name} subcommand expects ASQG input."
        )

    if first_format == "FASTG":
        _raise_logged_value_error(
            f"{file_name} looks like a FASTG file. "
            f"The {command_name} subcommand expects ASQG input."
        )

    _raise_logged_value_error(
        f"No ASQG segments were found in {file_name}. "
        f"The {command_name} subcommand expects an ASQG file containing VT records."
    )
