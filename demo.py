from agtools.graph_utils import get_graph

def main():
    assembly_graph_file = "tests/data/assembly_graph_with_scaffolds.gfa"
    contig_paths_file = "tests/data/contigs.paths"

    (
        assembly_graph, 
        contigs_map, 
        contig_names) = get_graph(
            assembly_graph_file, 
            contig_paths_file)
    
    print(f"Total number of vertices in the assembly graph: {len(list(assembly_graph.vs))}")
    print(f"Total number of edges in the assembly graph: {len(list(assembly_graph.es))}")

if __name__ == "__main__":
    main()