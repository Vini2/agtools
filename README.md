# agtools: Tools for manipulating assembly graphs

This repo contains useful scripts to manipulate assembly graph files from different assemblers.

Currently, I have added code for [SPAdes](https://github.com/ablab/spades) (including metaSPAdes) assembly graphs.

## Requirements

You should have Python and the following packages installed.

* [collections](https://docs.python.org/3/library/collections.html) - Usually installed by default with Python
* [re](https://docs.python.org/3/library/re.html) - Usually installed by default with Python
* [python-igraph](https://python.igraph.org/en/stable/index.html)

## Demo

Please see `demo.py` for an example.

```python
from agtools.spades_graph_utils import get_graph

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
        # Get neighbours of each contig
        neighbours = assembly_graph.neighbors(i, mode="all")

        # Print ID, contig label and number of neighbours
        print(
            assembly_graph.vs[i]["id"], assembly_graph.vs[i]["label"], len(neighbours)
        )


if __name__ == "__main__":
    main()
```

You will get `assembly_graph` which is an `igraph.Graph` object. You can get more details from the official [`python-igraph` documentation](https://python.igraph.org/en/stable/analysis.html#graph-analysis).