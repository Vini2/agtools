#!/usr/bin/env python3

import pytest

from agtools.assemblers import flye, megahit, myloasm


def test_flye_segment_path_parser_handles_path_normalization(tmp_path):
    contig_paths = tmp_path / "assembly_info.txt"
    contig_paths.write_text(
        "# comment\n"
        "seq_name\tsome_header\n"
        "contigA\t*,1,-2,*\n"
    )

    segment_contigs, contig_names, contig_name_to_id = (
        flye._get_segment_paths_and_contig_mapping(
            str(contig_paths), {"edge_1": 0, "edge_2": 1}
        )
    )

    assert contig_names == ["contigA"]
    assert contig_name_to_id == {"contigA": 0}
    assert segment_contigs[0] == {0}
    assert segment_contigs[1] == {0}


def test_myloasm_link_parser_filters_segments_not_in_fasta_index(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text(
        "S\tkeep\tAAAA\n"
        "S\tskip\tTTTT\n"
        "L\tkeep\t+\tskip\t+\t1M\n"
        "L\tkeep\t+\tkeep\t+\t1M\n"
    )

    contig_names, contig_name_to_id, edge_list, self_loops, lcount = (
        myloasm._get_links_and_contig_mapping_myloasm(str(gfa_file), {"keep": 0})
    )

    assert contig_names == ["keep"]
    assert contig_name_to_id == {"keep": 0}
    assert edge_list == []
    assert self_loops == [0]
    assert lcount == 2


def test_megahit_contig_graph_reports_missing_length_mapping(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tg1\tAAA\nS\tg2\tCCCC\n")

    contigs_file = tmp_path / "contigs.fa"
    contigs_file.write_text(">c1\nAAA\n>c2\nTTTTT\n")

    contig_graph = megahit.get_contig_graph(str(gfa_file), str(contigs_file))

    assert contig_graph.graph_to_contig_map["g1"] == "c1"
    assert str(contig_graph.get_contig_sequence("g1")) == "AAA"
    with pytest.raises(KeyError):
        contig_graph.get_contig_sequence("g2")

