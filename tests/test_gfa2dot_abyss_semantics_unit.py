#!/usr/bin/env python3

from agtools.commands.gfa2dot import gfa2dot


def test_gfa2dot_abyss_follows_oriented_vertex_and_edge_semantics(tmp_path):
    gfa_file = tmp_path / "graph.gfa"
    gfa_file.write_text("S\tA\tATGC\nS\tB\tGGGG\nL\tA\t+\tB\t-\t3M\n")

    target = tmp_path / "graph.gv"
    output_file = gfa2dot(str(gfa_file), abyss=True, output_path=str(target))
    content = target.read_text()

    assert output_file == str(target)

    # Each sequence has twin oriented vertices with length property.
    assert '"A+" [l=4]' in content
    assert '"A-" [l=4]' in content
    assert '"B+" [l=4]' in content
    assert '"B-" [l=4]' in content

    # Overlap edge and its reverse-complement twin with negative distance.
    assert '"A+" -> "B-" [d=-3]' in content
    assert '"B+" -> "A-" [d=-3]' in content

    # ABySS DOT uses oriented sequence names, not numeric vertex IDs.
    assert '"0" ->' not in content
    assert '"1" ->' not in content
