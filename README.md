# agtools: Tools for manipulating assembly graphs

This repo contains useful scripts to manipulate assembly graph files from different assemblers.

Currently, I have added code for [SPAdes](https://github.com/ablab/spades) (including metaSPAdes) assembly graphs.

Please see `demo.py` for an example.

```python
from agtools.spades_graph_utils import get_graph

def main():
    assembly_graph_file = "tests/data/assembly_graph_with_scaffolds.gfa"
    contig_paths_file = "tests/data/contigs.paths"

    (
        assembly_graph, 
        contigs_map, 
        contig_names) = get_graph(
            assembly_graph_file, 
            contig_paths_file)
    
    print(f"Total number of vertices in the assembly graph: {assembly_graph.vcount()}")
    print(f"Total number of edges in the assembly graph: {assembly_graph.ecount()}")

if __name__ == "__main__":
    main()
```

You will get `assembly_graph` which is an `igraph.Graph` object. You can get more details from the official [`python-igraph` documentation](https://python.igraph.org/en/stable/analysis.html#graph-analysis).