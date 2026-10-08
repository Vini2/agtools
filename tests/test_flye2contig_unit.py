#!/usr/bin/env python3

import pytest
from click.testing import CliRunner

from agtools.assemblers.flye import get_contig_links, get_contig_paths
from agtools.cli import main
from agtools.commands.flye2contig import _read_assembly_info, flye2contig

# A toy Flye-style edge graph.
#
#   contig_1 walks edge_1+ then edge_2-            (two edges merged into one contig)
#   contig_2 walks edge_3+                         (circular: edge_3 links to itself)
#   contig_3 walks edge_4+                         (joined to the end of contig_1)
#   contig_4 walks edge_5+                         (dropped: absent from the FASTA)
EXAMPLE_GFA = (
    "H\tVN:Z:1.0\n"
    "S\tedge_1\tAAAACCCC\tdp:i:10\n"
    "S\tedge_2\tGGGGTTTT\tdp:i:10\n"
    "S\tedge_3\tACGTACGT\tdp:i:30\n"
    "S\tedge_4\tTTTTAAAA\tdp:i:12\n"
    "S\tedge_5\tCCCCGGGG\tdp:i:3\n"
    "L\tedge_1\t+\tedge_2\t-\t0M\n"
    "L\tedge_2\t-\tedge_4\t+\t0M\n"
    "L\tedge_3\t+\tedge_3\t+\t0M\n"
    "P\tcontig_1\tedge_1+,edge_2-\t*\n"
    "P\tcontig_2\tedge_3+\t*\n"
    "P\tcontig_3\tedge_4+\t*\n"
    "P\tcontig_4\tedge_5+\t*\n"
)

EXAMPLE_FASTA = (
    ">contig_1\nAAAACCCCAAAACCCC\n>contig_2\nACGTACGT\n>contig_3\nTTTTAAAA\n"
)

EXAMPLE_INFO = (
    "#seq_name\tlength\tcov.\tcirc.\trepeat\tmult.\talt_group\tgraph_path\n"
    "contig_1\t16\t10\tN\tN\t1\t*\t1,-2\n"
    "contig_2\t8\t30\tY\tN\t1\t*\t3\n"
    "contig_3\t8\t12\tN\tN\t1\t*\t4\n"
)


@pytest.fixture
def flye_assembly(tmp_path):
    gfa_file = tmp_path / "assembly_graph.gfa"
    gfa_file.write_text(EXAMPLE_GFA)

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file.write_text(EXAMPLE_FASTA)

    info_file = tmp_path / "assembly_info.txt"
    info_file.write_text(EXAMPLE_INFO)

    return gfa_file, fasta_file, info_file


def _records(content, tag):
    return [
        line.split("\t") for line in content.splitlines() if line.split("\t")[0] == tag
    ]


def test_get_contig_paths_parses_oriented_steps(flye_assembly):
    gfa_file, _, _ = flye_assembly

    paths = get_contig_paths(str(gfa_file))

    assert paths["contig_1"] == [("edge_1", "+"), ("edge_2", "-")]
    assert paths["contig_2"] == [("edge_3", "+")]
    assert set(paths) == {"contig_1", "contig_2", "contig_3", "contig_4"}


def test_get_contig_paths_ignores_wildcard_steps(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tedge_1\tAC\nP\tcontig_1\t*,edge_1+,*\t*\n")

    assert get_contig_paths(str(gfa_file)) == {"contig_1": [("edge_1", "+")]}


def test_get_contig_links_joins_outward_tips_only(flye_assembly):
    gfa_file, _, _ = flye_assembly

    paths = get_contig_paths(str(gfa_file))
    links = get_contig_links(str(gfa_file), paths)

    # edge_1/edge_2 are joined inside contig_1, so that link is not lifted.
    assert ("contig_1", "+", "contig_3", "+", "0M") in links
    assert ("contig_2", "+", "contig_2", "+", "0M") in links
    assert len(links) == 2


def test_get_contig_links_respects_the_contig_subset(flye_assembly):
    gfa_file, _, _ = flye_assembly

    paths = get_contig_paths(str(gfa_file))
    paths.pop("contig_3")

    links = get_contig_links(str(gfa_file), paths)

    assert all("contig_3" not in link for link in links)


def test_read_assembly_info_extracts_length_coverage_and_topology(flye_assembly):
    _, _, info_file = flye_assembly

    info = _read_assembly_info(str(info_file))

    assert info["contig_2"] == {"length": 8, "coverage": 30, "circular": True}
    assert info["contig_1"]["circular"] is False


def test_flye2contig_writes_contig_level_graph(flye_assembly, tmp_path):
    gfa_file, fasta_file, info_file = flye_assembly
    output = tmp_path / "contigs.gfa"

    output_file = flye2contig(
        str(gfa_file), str(fasta_file), str(info_file), str(output)
    )
    content = output.read_text()

    assert output_file == str(output)
    assert content.startswith("H\tVN:Z:1.0\n")

    segments = _records(content, "S")
    assert [record[1] for record in segments] == ["contig_1", "contig_2", "contig_3"]
    assert segments[0][2] == "AAAACCCCAAAACCCC"
    assert "LN:i:16" in segments[0]
    assert "dp:i:10" in segments[0]
    assert "TP:Z:linear" in segments[0]
    assert "TP:Z:circular" in segments[1]

    links = {tuple(record[1:6]) for record in _records(content, "L")}
    assert links == {
        ("contig_1", "+", "contig_3", "+", "0M"),
        ("contig_2", "+", "contig_2", "+", "0M"),
    }


def test_flye2contig_drops_paths_missing_from_the_fasta(flye_assembly, tmp_path):
    gfa_file, fasta_file, info_file = flye_assembly
    output = tmp_path / "contigs.gfa"

    flye2contig(str(gfa_file), str(fasta_file), str(info_file), str(output))

    assert "contig_4" not in output.read_text()
    assert "edge_" not in output.read_text()


def test_flye2contig_works_without_assembly_info(flye_assembly, tmp_path):
    gfa_file, fasta_file, _ = flye_assembly
    output = tmp_path / "contigs.gfa"

    flye2contig(str(gfa_file), str(fasta_file), None, str(output))
    content = output.read_text()

    assert "LN:i:16" in content
    assert "dp:i:" not in content
    assert "TP:Z:" not in content


def test_flye2contig_rejects_graphs_without_contig_paths(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("H\tVN:Z:1.0\nS\tedge_1\tACGT\n")

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file.write_text(">contig_1\nACGT\n")

    with pytest.raises(ValueError, match="No contig paths"):
        flye2contig(
            str(gfa_file),
            str(fasta_file),
            None,
            str(tmp_path / "contigs.gfa"),
        )


def test_flye2contig_rejects_non_gfa_input(tmp_path):
    fastg_file = tmp_path / "graph.fastg"
    fastg_file.write_text(">NODE_1_length_4\nACGT\n")

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file.write_text(">contig_1\nACGT\n")

    with pytest.raises(ValueError, match="FASTG"):
        flye2contig(
            str(fastg_file),
            str(fasta_file),
            None,
            str(tmp_path / "contigs.gfa"),
        )


def test_flye2contig_cli_writes_output(flye_assembly, tmp_path):
    gfa_file, fasta_file, info_file = flye_assembly
    output = tmp_path / "contigs.gfa"

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "flye2contig",
            "-g",
            str(gfa_file),
            "-f",
            str(fasta_file),
            "-i",
            str(info_file),
            "-o",
            str(output),
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0
    assert output.exists()
    assert "contig_1" in output.read_text()


# ---------------------------------------------------------------------------
# Orientation
#
# A GFA link ``L u ou v ov`` joins the exit of ``u`` traversed in ``ou`` to the
# entry of ``v`` traversed in ``ov``. Lifting it to contigs must preserve that:
# walking contig A in orientation ``oa`` has to leave through ``u^ou``, and
# walking contig B in ``ob`` has to arrive through ``v^ov``. The cases below
# are derived by hand from that rule.
# ---------------------------------------------------------------------------

ORIENTATION_GFA = (
    "H\tVN:Z:1.0\n"
    "S\tedge_1\tAAAACCCC\n"
    "S\tedge_2\tGGGGTTTT\n"
    "S\tedge_3\tACGTACGT\n"
    "L\tedge_1\t+\tedge_2\t-\t0M\n"
    "L\tedge_2\t+\tedge_3\t+\t0M\n"
)

ORIENTATION_FASTA = ">contig_1\nAAAACCCC\n>contig_2\nGGGGTTTT\n>contig_3\nACGTACGT\n"


def _run_orientation_case(tmp_path, path_records):
    gfa_file = tmp_path / "assembly_graph.gfa"
    gfa_file.write_text(ORIENTATION_GFA + path_records)

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file.write_text(ORIENTATION_FASTA)

    output = tmp_path / "contigs.gfa"
    flye2contig(str(gfa_file), str(fasta_file), None, str(output))

    return {tuple(record[1:5]) for record in _records(output.read_text(), "L")}


def test_links_follow_the_orientation_of_the_edge_link(tmp_path):
    # contig_1 exits through edge_1+, contig_2 is entered through edge_2-.
    # "L edge_1 + edge_2 -" therefore lifts to "L contig_1 + contig_2 +",
    # because walking contig_2 forwards enters through edge_2-.
    links = _run_orientation_case(
        tmp_path,
        "P\tcontig_1\tedge_1+\t*\nP\tcontig_2\tedge_2-\t*\n",
    )

    assert links == {("contig_1", "+", "contig_2", "+")}


def test_flipping_a_path_orientation_flips_the_link_sign(tmp_path):
    # Same graph, but contig_2 now walks edge_2 forwards. Entering through
    # edge_2- now means traversing contig_2 in reverse, so the link sign flips.
    links = _run_orientation_case(
        tmp_path,
        "P\tcontig_1\tedge_1+\t*\nP\tcontig_2\tedge_2+\t*\n",
    )

    assert links == {("contig_1", "+", "contig_2", "-")}


def test_flipping_the_source_path_flips_the_source_sign(tmp_path):
    # contig_1 now walks edge_1 in reverse, so leaving through edge_1+ means
    # traversing contig_1 backwards.
    links = _run_orientation_case(
        tmp_path,
        "P\tcontig_1\tedge_1-\t*\nP\tcontig_2\tedge_2-\t*\n",
    )

    assert links == {("contig_1", "-", "contig_2", "+")}


def test_links_attach_to_path_ends_not_to_interior_edges(tmp_path):
    # contig_2 walks edge_2+ then edge_3+.
    #
    #   "L edge_2 + edge_3 +" is exactly the junction interior to contig_2
    #   (it joins the 3' end of edge_2 to the 5' end of edge_3, both buried
    #   inside the contig), so it must NOT be lifted.
    #
    #   "L edge_1 + edge_2 +" lands on the 5' end of edge_2, which is the
    #   outward-facing start of contig_2, so it must be lifted. Leaving
    #   contig_2 through its start means traversing it in reverse, so that tip
    #   belongs to contig_2-, and entering it there is contig_2+.
    gfa_file = tmp_path / "assembly_graph.gfa"
    gfa_file.write_text(
        "H\tVN:Z:1.0\n"
        "S\tedge_1\tAAAACCCC\n"
        "S\tedge_2\tGGGGTTTT\n"
        "S\tedge_3\tACGTACGT\n"
        "L\tedge_1\t+\tedge_2\t+\t0M\n"
        "L\tedge_2\t+\tedge_3\t+\t0M\n"
        "P\tcontig_1\tedge_1+\t*\n"
        "P\tcontig_2\tedge_2+,edge_3+\t*\n"
    )

    fasta_file = tmp_path / "assembly.fasta"
    fasta_file.write_text(">contig_1\nAAAACCCC\n>contig_2\nGGGGTTTTACGTACGT\n")

    output = tmp_path / "contigs.gfa"
    flye2contig(str(gfa_file), str(fasta_file), None, str(output))

    links = {tuple(record[1:5]) for record in _records(output.read_text(), "L")}

    assert links == {("contig_1", "+", "contig_2", "+")}


def test_every_link_is_backed_by_an_orientation_matching_edge_link(
    flye_assembly, tmp_path
):
    """Re-derive each link from the P records, without the tip helper."""

    gfa_file, fasta_file, info_file = flye_assembly
    output = tmp_path / "contigs.gfa"

    flye2contig(str(gfa_file), str(fasta_file), str(info_file), str(output))

    paths = get_contig_paths(str(gfa_file))
    edge_links = {
        (record[1], record[2], record[3], record[4])
        for record in _records(EXAMPLE_GFA, "L")
    }

    def flip(orientation):
        return "-" if orientation == "+" else "+"

    def walk(contig, orientation):
        steps = paths[contig]
        if orientation == "+":
            return steps
        return [(edge, flip(step)) for edge, step in reversed(steps)]

    for record in _records(output.read_text(), "L"):
        from_contig, from_orient, to_contig, to_orient = record[1:5]
        exit_edge, exit_orient = walk(from_contig, from_orient)[-1]
        entry_edge, entry_orient = walk(to_contig, to_orient)[0]

        forward = (exit_edge, exit_orient, entry_edge, entry_orient)
        reverse = (entry_edge, flip(entry_orient), exit_edge, flip(exit_orient))

        assert forward in edge_links or reverse in edge_links, record
