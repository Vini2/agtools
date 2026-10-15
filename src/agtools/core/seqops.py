#!/usr/bin/env python3
"""Lightweight nucleotide sequence operations.

These replace the handful of ``Bio.Seq`` helpers ``agtools`` used to depend
on, so that sequences can be passed around as plain :class:`str` objects.
"""

__all__ = ["complement", "reverse_complement"]

# Full IUPAC nucleotide complement table, preserving case and passing any
# other character (e.g. gap characters) through unchanged.
_COMPLEMENT = str.maketrans(
    "ACGTUMRWSYKVHDBNacgtumrwsykvhdbn",
    "TGCAAKYWSRMBDHVNtgcaakywsrmbdhvn",
)


def complement(sequence: str) -> str:
    """
    Return the complement of a nucleotide sequence.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence. IUPAC ambiguity codes are supported and case is
        preserved; unrecognised characters pass through unchanged.

    Returns
    -------
    str
        The complement of the input sequence.

    Examples
    --------
    >>> complement("ATGC")
    'TACG'
    """
    return str(sequence).translate(_COMPLEMENT)


def reverse_complement(sequence: str) -> str:
    """
    Return the reverse complement of a nucleotide sequence.

    Parameters
    ----------
    sequence : str
        Nucleotide sequence. IUPAC ambiguity codes are supported and case is
        preserved; unrecognised characters pass through unchanged.

    Returns
    -------
    str
        The reverse complement of the input sequence.

    Examples
    --------
    >>> reverse_complement("ATGC")
    'GCAT'
    """
    return str(sequence).translate(_COMPLEMENT)[::-1]
