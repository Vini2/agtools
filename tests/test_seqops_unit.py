#!/usr/bin/env python3

from agtools.core.seqops import complement, reverse_complement


def test_complement_handles_case_and_ambiguity_codes():
    assert complement("ATGC") == "TACG"
    assert complement("atgc") == "tacg"
    assert complement("RYKMSWN") == "YRMKSWN"


def test_reverse_complement_matches_expected():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AAAACCCGGT") == "ACCGGGTTTT"
    assert reverse_complement("") == ""


def test_reverse_complement_is_an_involution():
    sequence = "ACGTTGCANNRY"
    assert reverse_complement(reverse_complement(sequence)) == sequence


def test_unknown_characters_pass_through():
    assert reverse_complement("AT-GC") == "GC-AT"
