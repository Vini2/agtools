#!/usr/bin/env python3

import pytest

from agtools.commands.concat import concat


def test_concat_orders_known_tags_and_appends_unknown_tags(tmp_path):
    graph_1 = tmp_path / "graph1.gfa"
    graph_2 = tmp_path / "graph2.gfa"

    graph_1.write_text(
        "W\twalk1\t*\t*\t*\t>seg1\n"
        "S\tseg1\tAT\n"
        "H\tVN:Z:1.0\n"
        "X\tcustom1\n"
    )
    graph_2.write_text(
        "P\tpath1\tseg1+\t*\n" "S\tseg2\tGC\n" "L\tseg1\t+\tseg2\t+\t1M\n" "#\tcomment\n"
    )

    target = tmp_path / "concatenated_graph.gfa"
    output_file = concat([str(graph_1), str(graph_2)], str(target))
    lines = target.read_text().splitlines()

    assert output_file == str(tmp_path / "concatenated_graph.gfa")
    assert lines.index("#\tcomment") < lines.index("H\tVN:Z:1.0")
    assert lines.index("H\tVN:Z:1.0") < lines.index("S\tseg1\tAT")
    assert lines.index("S\tseg1\tAT") < lines.index("S\tseg2\tGC")
    assert lines.index("S\tseg2\tGC") < lines.index("L\tseg1\t+\tseg2\t+\t1M")
    assert lines.index("L\tseg1\t+\tseg2\t+\t1M") < lines.index("P\tpath1\tseg1+\t*")
    assert lines.index("P\tpath1\tseg1+\t*") < lines.index("W\twalk1\t*\t*\t*\t>seg1")
    assert lines[-1] == "X\tcustom1"


@pytest.mark.parametrize(
    ("line_1", "line_2"),
    [
        ("S\tseg1\tATGC\n", "S\tseg1\tGCTA\n"),
        ("P\tpath1\tseg1+\t*\n", "P\tpath1\tseg2+\t*\n"),
        ("W\twalk1\t*\t*\t*\t>seg1\n", "W\twalk1\t*\t*\t*\t>seg2\n"),
    ],
)
def test_concat_raises_on_duplicate_ids(tmp_path, line_1, line_2):
    graph_1 = tmp_path / "graph1.gfa"
    graph_2 = tmp_path / "graph2.gfa"
    graph_1.write_text(line_1)
    graph_2.write_text(line_2)

    target = tmp_path / "concatenated_graph.gfa"
    with pytest.raises(SystemExit) as exc:
        concat([str(graph_1), str(graph_2)], str(target))

    assert exc.value.code == 1
