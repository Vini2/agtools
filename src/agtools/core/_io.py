#!/usr/bin/env python3
"""Transparent compressed/plain text I/O helpers.

All ``agtools`` readers of line-oriented text formats (GFA, FASTG, ASQG,
assembler path files) go through :func:`open_text` so that both plain and
gzip-compressed inputs are accepted without the caller caring which it has.

FASTA/FASTQ reading and writing is *not* handled here -- that is delegated to
`rsbio-seq <https://github.com/anuradhawick/rsbio-seq>`_, which has a Rust
core and handles plain, gzip and BGZF input natively.
"""

import gzip
import io

__all__ = [
    "is_gzip",
    "is_bgzf",
    "open_text",
    "strip_compression_suffix",
]

_GZIP_MAGIC = b"\x1f\x8b"
_FEXTRA = 0x04


def is_gzip(file_path) -> bool:
    """
    Check whether a file is gzip-compressed.

    Detection is based on the file's magic bytes rather than its name, so
    gzip files without a ``.gz`` suffix are still recognised.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the file to inspect.

    Returns
    -------
    bool
        True if the file starts with the gzip magic bytes.

    Examples
    --------
    >>> is_gzip("assembly_graph.gfa.gz")
    True
    """
    try:
        with open(file_path, "rb") as handle:
            return handle.read(2) == _GZIP_MAGIC
    except OSError:
        return False


def is_bgzf(file_path) -> bool:
    """
    Check whether a file is BGZF-compressed (``bgzip`` output).

    BGZF is a gzip-compatible block format. Unlike ordinary gzip it supports
    random access, which is what allows ``rsbio-seq`` to serve indexed
    (``.fai`` + ``.gzi``) lookups from a compressed FASTA file.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the file to inspect.

    Returns
    -------
    bool
        True if the file is a BGZF stream, False for plain files and for
        ordinary gzip files.

    Examples
    --------
    >>> is_bgzf("contigs.fa.gz")  # produced with `bgzip contigs.fa`
    True
    """
    try:
        with open(file_path, "rb") as handle:
            header = handle.read(12)
            if len(header) < 12 or header[:2] != _GZIP_MAGIC:
                return False
            if not header[3] & _FEXTRA:
                return False
            extra_len = int.from_bytes(header[10:12], "little")
            extra = handle.read(extra_len)
    except OSError:
        return False

    pos = 0
    while pos + 4 <= len(extra):
        sub_id = extra[pos : pos + 2]
        sub_len = int.from_bytes(extra[pos + 2 : pos + 4], "little")
        if sub_id == b"BC":
            return True
        pos += 4 + sub_len

    return False


def open_text(file_path, mode: str = "rt", buffering: int = io.DEFAULT_BUFFER_SIZE):
    """
    Open a text file for reading or writing, transparently handling gzip.

    For reading, compression is auto-detected from the file's magic bytes.
    For writing, compression is chosen from a ``.gz`` suffix on the path.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to the file to open.
    mode : str, optional
        Text mode, default ``"rt"``. Binary modes are not supported.
    buffering : int, optional
        Buffer size for plain-text files. Ignored for gzip streams, which
        manage their own buffering.

    Returns
    -------
    io.TextIOBase
        An open text handle. The handle is seekable for both plain and gzip
        inputs, so byte-offset based random access keeps working -- though
        seeking into an ordinary gzip stream is much slower, since the data
        before the target offset has to be decompressed again.

    Examples
    --------
    >>> with open_text("assembly_graph.gfa.gz") as handle:
    ...     first_line = handle.readline()
    """
    if "b" in mode:
        raise ValueError("open_text only supports text modes")

    writing = any(flag in mode for flag in ("w", "a", "x"))
    compressed = str(file_path).endswith(".gz") if writing else is_gzip(file_path)

    if compressed:
        return gzip.open(file_path, mode)

    return io.open(file_path, mode=mode, buffering=buffering)


def strip_compression_suffix(file_path) -> str:
    """
    Remove a trailing ``.gz`` from a path, if present.

    Useful when deriving an output name from a compressed input name.

    Parameters
    ----------
    file_path : str or pathlib.Path
        Path to strip.

    Returns
    -------
    str
        The path without its ``.gz`` suffix.

    Examples
    --------
    >>> strip_compression_suffix("contigs.fasta.gz")
    'contigs.fasta'
    """
    text = str(file_path)
    return text[:-3] if text.endswith(".gz") else text
