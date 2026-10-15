#!/usr/bin/env python3

import struct
import zlib

import pytest

# Standard 28-byte BGZF end-of-file marker block.
BGZF_EOF = (
    b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00"
    b"\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
)


def _write_bgzf(path, text: str):
    """Write ``text`` as a single-block BGZF file, as ``bgzip`` would."""
    payload = text.encode()
    compressor = zlib.compressobj(9, zlib.DEFLATED, -15)
    compressed = compressor.compress(payload) + compressor.flush()
    block_size = len(compressed) + 26

    header = (
        b"\x1f\x8b\x08\x04"
        + b"\x00" * 6
        + struct.pack("<H", 6)
        + b"BC"
        + struct.pack("<HH", 2, block_size - 1)
    )
    trailer = struct.pack("<II", zlib.crc32(payload) & 0xFFFFFFFF, len(payload))

    path.write_bytes(header + compressed + trailer + BGZF_EOF)
    return path


@pytest.fixture
def write_bgzf():
    """Return a helper that writes BGZF-compressed text to a given path."""
    return _write_bgzf
