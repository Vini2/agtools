# API Tutorial

## Loading graphs from the Python package interface

### Loading a GFA file

```python
from agtools.core.unitig_graph import UnitigGraph

graph_file = "tests/data/ESC/assembly_graph_with_scaffolds.gfa"
unitig_graph = UnitigGraph.from_gfa(graph_file)

print(f"Graph loaded from file: {unitig_graph.file_path}")
print(f"Number of segments: {unitig_graph.vcount}")
print(f"Number of links: {unitig_graph.ecount}")
```

### Loading a SPAdes graph

```python
from agtools.assemblers import spades

graph_file = "tests/data/ESC/assembly_graph_with_scaffolds.gfa"
contigs_file = "tests/data/ESC/contigs.fasta"
contig_paths_file = "tests/data/ESC/contigs.paths"

contig_graph = spades.get_contig_graph(graph_file, contigs_file, contig_paths_file)
unitig_graph = spades.get_unitig_graph(graph_file)
```

### Loading a MEGAHIT graph

```python
from agtools.assemblers import megahit

graph_file = "tests/data/5G/final.gfa"
contig_file = "tests/data/5G/final.contigs.fa"

contig_graph = megahit.get_contig_graph(graph_file, contig_file)
```

### Loading a Flye graph

```python
from agtools.assemblers import flye

graph_file = "tests/data/1Y3B/assembly_graph.gfa"
contig_paths_file = "tests/data/1Y3B/assembly_info.txt"
contigs_file = "tests/data/1Y3B/assembly.fasta"

contig_graph = flye.get_contig_graph(graph_file, contigs_file, contig_paths_file)
unitig_graph = flye.get_unitig_graph(graph_file)
```