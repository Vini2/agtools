## Obtaining stats of an assembly graph

*agtools* can generate basic graph and sequence-based statistics for an assembly graph. You can use the `stats` subcommand provided through the command-line interface. Please refer to the [CLI reference](../cli.md) for further details on the `stats` subcommand.

```bash
agtools stats -g assembly_graph_with_scaffolds.gfa -o results/graph_stats.txt
```

A text file with the following details will be generated.

```text
Basic graph statistics for tests/data/ESC/assembly_graph_with_scaffolds.gfa:
Number of segments: 982
Number of links: 1265
Number of self-loops: 1
Number of connected components: 4
Average node degree: 2

Sequence-based statistics for tests/data/ESC/assembly_graph_with_scaffolds.gfa:
Total length of segments: 8337494 bp
Average segment length: 8490 bp
N50: 60706 bp
L50: 36 segment(s)
GC content: 43.52%
```
