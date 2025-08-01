# Command-Line Interface (CLI) Reference

agtools provides tools for manipulating assembly graphs.

Run `agtools --help` or `agtools -h` to list the help message for agtools.

```bash
Usage: agtools [OPTIONS] COMMAND [ARGS]...

  agtools: Tools for manipulating assembly graphs

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  stats      Compute statistics about the graph
  rename     Rename segments, paths and walks in a GFA file
  concat     Concatenate two or more GFA files
  filter     Filter segments from GFA file
  clean      Clean a GFA file based on segments in a FASTA file
  component  Extract a component containing a given segment
  fastg2gfa  Convert FASTG file to GFA format
  asqg2gfa   Convert ASQG file to GFA format
  gfa2fastg  Convert GFA file to FASTG format
  gfa2dot    Convert GFA file to DOT format (GraphViz)
  gfa2fasta  Get segments in FASTA format
  gfa2adj    Get adjacency matrix of the assembly graph
```

## `stats`

Compute statistics about the graph.

Run `agtools stats --help` or `agtools stats -h` to list the help message for agtools stats.

```bash
Usage: agtools stats [OPTIONS]

  Compute statistics about the graph

Options:
  -g, --graph PATH   path(s) to the assembly graph file(s)  [required]
  -o, --output PATH  path to the output folder  [r
```