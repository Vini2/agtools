# agtools Command-Line Interface Reference

*agtools* provides tools for manipulating assembly graphs.

Run `agtools --help` or `agtools -h` to list the help message for *agtools*.

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
  -o, --output PATH  path to the output folder  [required]
  -h, --help         Show this message and exit.
```

**Inputs**

* Assembly graph file in GFA format

**Outputs**

* Text file with the statistics

## `rename`

Rename segments, paths and walks in a GFA file.

Run `agtools rename --help` or `agtools rename -h` to list the help message for agtools rename.

```bash
Usage: agtools rename [OPTIONS]

  Rename segments, paths and walks in a GFA file

Options:
  -g, --graph PATH   path(s) to the assembly graph file(s)  [required]
  -p, --prefix TEXT  prefix for the graph elements  [default: ""]
  -o, --output PATH  path to the output folder  [required]
  -h, --help         Show this message and exit.
```

**Inputs**

* Assembly graph file in GFA format
* Prefix to prepend segments, paths and walks in a GFA file

**Outputs**

* Assembly graph file in GFA format with renamed elements

## `concat`

Concatenate two or more GFA files.

Run `agtools concat --help` or `agtools concat -h` to list the help message for agtools concat.

```bash
Usage: agtools concat [OPTIONS]

  Concatenate two or more GFA files

Options:
  -g, --graph PATH   path(s) to the assembly graph file(s)  [required]
  -o, --output PATH  path to the output folder  [required]
  -h, --help         Show this message and exit.
```

**Inputs**

* Assembly graph files in GFA format. You can provide each file followed by `-g`

**Outputs**

* Concatenated assembly graph file in GFA format

## `filter`

Filter segments from a GFA file based on length.

Run `agtools filter --help` or `agtools filter -h` to list the help message for agtools filter.

```bash
Usage: agtools filter [OPTIONS]

  Filter segments from GFA file

Options:
  -g, --graph PATH          path(s) to the assembly graph file(s)  [required]
  -l, --min-length INTEGER  minimum length of segments to keep  [default: 100;
                            required]
  -o, --output PATH         path to the output folder  [required]
  -h, --help                Show this message and exit.
```

**Inputs**

* Assembly graph file in GFA format
* Minimum length of segments to keep

**Outputs**

* Filtered assembly graph file in GFA format

## `clean`

Clean a GFA file based on segments in a FASTA file.

Run `agtools clean --help` or `agtools clean -h` to list the help message for agtools clean.

```bash
Usage: agtools clean [OPTIONS]

  Clean a GFA file based on segments in a FASTA file

Options:
  -g, --graph PATH      path(s) to the assembly graph file(s)  [required]
  -f, --fasta PATH      path to the FASTA file
  -a, --assembler TEXT  assembler name (if assembler used is myloasm)
  -o, --output PATH     path to the output folder  [required]
  -h, --help            Show this message and exit.
```

**Inputs**

* Assembly graph file in GFA format
* FASTA file with segment sequences
* Name of the assembler

**Outputs**

* Cleaned assembly graph file in GFA format. Removes segments which are not found in the FASTA file and adds the sequences if they are not present in the GFA file (to be compatible with the [GFA specification](https://gfa-spec.github.io/GFA-spec/GFA1.html)).

## `component`

Extract a component containing a given segment.

Run `agtools component --help` or `agtools component -h` to list the help message for agtools component.

```bash
Usage: agtools component [OPTIONS]

  Extract a component containing a given segment

Options:
  -g, --graph PATH    path(s) to the assembly graph file(s)  [required]
  -s, --segment TEXT  segment ID  [required]
  -o, --output PATH   path to the output folder  [required]
  -h, --help          Show this message and exit.
```

**Inputs**

* Assembly graph file in GFA format
* Segment ID of interest

**Outputs**

* Assembly graph file of the constituent component in GFA format

## `fastg2gfa`

Convert a FASTG file to GFA format.

Run `agtools fastg2gfa --help` or `agtools fastg2gfa -h` to list the help message for agtools fastg2gfa.

```bash
Usage: agtools fastg2gfa [OPTIONS]

  Convert FASTG file to GFA format

Options:
  -g, --graph PATH     path(s) to the assembly graph file(s)  [required]
  -k, --ksize INTEGER  k-mer size used for the assembly  [default: 141;
                       required]
  -o, --output PATH    path to the output folder  [required]
  -h, --help           Show this message and exit.
```

**Inputs**

* Assembly graph file in FASTG format
* k-mer size used for the assembly

**Outputs**

* Assembly graph file in GFA format

## `asqg2gfa`

Convert a ASQG file to GFA format.

Run `agtools asqg2gfa --help` or `agtools asqg2gfa -h` to list the help message for agtools asqg2gfa.

```bash
Usage: agtools asqg2gfa [OPTIONS]

  Convert ASQG file to GFA format

Options:
  -g, --graph PATH   path(s) to the assembly graph file(s)  [required]
  -o, --output PATH  path to the output folder  [required]
  -h, --help         Show this message and exit.
```

**Inputs**

* Assembly graph file in ASQG format

**Outputs**

* Assembly graph file in GFA format

## `gfa2dot`

Convert GFA file to GraphViz DOT format or ABySS DOT format.

Run `agtools gfa2dot --help` or `agtools gfa2dot -h` to list the help message for agtools gfa2dot.

```bash
Usage: agtools gfa2dot [OPTIONS]

  Convert GFA file to DOT format (GraphViz)

Options:
  -g, --graph PATH   path(s) to the assembly graph file(s)  [required]
  -ab, --abyss       use the ABySS DOT format for the output
  -o, --output PATH  path to the output folder  [required]
  -h, --help         Show this message and exit.
```

**Inputs**

* Assembly graph file in GFA format

**Outputs**

* Assembly graph file in [GraphViz DOT](https://graphviz.org/doc/info/lang.html) format or [ABySS DOT](https://github.com/bcgsc/abyss/wiki/ABySS-File-Formats#dot) format.

## `gfa2fasta`

Get segments from a GFA file in FASTA format

Run `agtools gfa2fasta --help` or `agtools gfa2fasta -h` to list the help message for agtools gfa2fasta.

```bash
Usage: agtools gfa2fasta [OPTIONS]

  Get segments in FASTA format

Options:
  -g, --graph PATH   path(s) to the assembly graph file(s)  [required]
  -o, --output PATH  path to the output folder  [required]
  -h, --help         Show this message and exit.
```

**Inputs**

* Assembly graph file in GFA format

**Outputs**

* Segment sequences in FASTA format

## `gfa2adj`

Get adjacency matrix of the assembly graph

Run `agtools gfa2adj --help` or `agtools gfa2adj -h` to list the help message for agtools gfa2adj.

```bash
Usage: agtools gfa2adj [OPTIONS]

  Get adjacency matrix of the assembly graph

Options:
  -g, --graph PATH         path(s) to the assembly graph file(s)  [required]
  --delimiter [comma|tab]  delimiter for adjacency file. Supports a comma and
                           a tab.  [default: comma]
  -o, --output PATH        path to the output folder  [required]
  -h, --help               Show this message and exit.
```

**Inputs**

* Assembly graph file in GFA format

**Outputs**

* A delimited file of the adjancency matrix