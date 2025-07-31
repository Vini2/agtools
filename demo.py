from agtools.assemblers.spades import get_graph


def main():
    assembly_graph_file = "tests/data/assembly_graph_with_scaffolds.gfa"
    contig_paths_file = "tests/data/contigs.paths"

    (assembly_graph, contigs_map, contig_names) = get_graph(
        assembly_graph_file, contig_paths_file
    )

    # Total number of vertices
    print(f"Total number of vertices in the assembly graph: {assembly_graph.vcount()}")

    # Total number of edges
    print(f"Total number of edges in the assembly graph: {assembly_graph.ecount()}")

    # Iterate through the contigs
    for i in range(assembly_graph.vcount()):
        # Get neighbors of each contig
        neighbors = assembly_graph.neighbors(i, mode="all")

        # Print ID, contig label and number of neighbors
        print(
            assembly_graph.vs[i]["id"], assembly_graph.vs[i]["label"], len(neighbors)
        )


if __name__ == "__main__":
    main()
