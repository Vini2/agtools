#!/usr/bin/env python3

from agtools.core.graph import parse_fastg

def _extract_links(segments, fixed_overlap) -> list:

    links = []

    for from_node, neighbors in segments.items():
        for raw_to_node in neighbors:
            if not raw_to_node:
                continue

            from_ori = "+" if not from_node.endswith("'") else "-"
            to_ori = "+" if not raw_to_node.endswith("'") else "-"

            from_clean = from_node.strip("'")
            to_clean = raw_to_node.strip("'")

            links.append((from_clean, from_ori, to_clean, to_ori, f"{fixed_overlap}M"))

    return links


def _write_gfa(segments, links, output_path) -> str:

    output_file = f"{output_path}/converted_graph.gfa"
    with open(output_file, 'w') as f:
        for seg_id, seq in segments.items():
            f.write(f'S\t{seg_id}\t{seq}\n')
        for from_node, from_ori, to_node, to_ori, overlap in links:
            f.write(f'L\t{from_node}\t{from_ori}\t{to_node}\t{to_ori}\t{overlap}\n')

    return output_file


def fastg2gfa(fastg_path, k_overlap, gfa_path) -> str:

    segments, edges = parse_fastg(fastg_path)
    links = _extract_links(edges, fixed_overlap=k_overlap)
    output_file = _write_gfa(segments, links, gfa_path)

    return output_file