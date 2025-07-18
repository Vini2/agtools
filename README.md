# agtools: Tools for manipulating assembly graphs

![GitHub License](https://img.shields.io/github/license/Vini2/agtools)
[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](#)
[![CI](https://github.com/Vini2/agtools/actions/workflows/testing_python_app.yml/badge.svg)](https://github.com/Vini2/agtools/actions/workflows/testing_python_app.yml)
[![codecov](https://codecov.io/gh/Vini2/agtools/graph/badge.svg?token=nYzx0Pd0h6)](https://codecov.io/gh/Vini2/agtools)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

`agtools` is a toolkit for manipulating assembly graphs, with a focus on the [Graphical Fragment Assembly (GFA) format](https://github.com/GFA-spec/GFA-spec). It offers a command-line interface for tasks such as graph format conversion, segment filtering, and component extraction. Additionally, it provides a Python package interface that exposes assembler-specific functionality for advanced analysis and integration.

## Requirements

You should have Python and the following packages installed.

* [flit](https://flit.pypa.io/en/stable/)
* [click](https://click.palletsprojects.com/en/stable/)
* [loguru](https://loguru.readthedocs.io/en/stable/)
* [bidict](https://bidict.readthedocs.io/en/main/intro.html)
* [python-igraph](https://python.igraph.org/en/stable/index.html)
* [biopython](https://biopython.org/)
* [pandas](https://pandas.pydata.org/)
* [collections](https://docs.python.org/3/library/collections.html) - Usually installed by default with Python
* [re](https://docs.python.org/3/library/re.html) - Usually installed by default with Python

## Installing `agtools`

### For development

Please follow the steps below to install `agtools` using `flit` for development.

```bash
# clone repository
git clone https://github.com/Vini2/agtools.git

# move to gbintk directory
cd agtools

# create and activate conda env
conda env create -f environment.yml
conda activate agtools

# install using flit
flit install -s --python `which python`

# test installation
agtools --help
```

## Available subcommands in `agtools`

Run `agtools --help` or `agtools -h` to list the help message for `agtools`.

```bash
Usage: agtools [OPTIONS] COMMAND [ARGS]...

  agtools: Tools for manipulating assembly graphs

Options:
  -v, --version  Show the version and exit.
  -h, --help     Show this message and exit.

Commands:
  stats      Compute statistics about the graph
  rename     Rename segments in a GFA file
  merge      Merge two or more GFA files
  filter     Filter segments from GFA file
  component  Extract a component containing a given segment
  fastg2gfa  Convert FASTG file to GFA format
  asqg2gfa   Convert ASQG file to GFA format
  gfa2fastg  Convert GFA file to FASTG format
  gfa2fasta  Get segments in FASTA format
  gfa2adj    Get adjacency matrix of the assembly graph
```

## Loading graphs from the Python package interface

### Loading a GFA file

```python
from agtools.core.graph import UnitigGraph

unitig_graph = UnitigGraph.from_gfa(graph_file)

print(f"Number of segments: {unitig_graph.graph.vcount()}")
print(f"Number of links: {unitig_graph.graph.ecount()}")
```

### Loading a SPAdes graph

```python
from agtools.assemblers import spades

graph_file = "tests/data/ESC/assembly_graph_with_scaffolds.gfa"
contig_paths_file = "tests/data/ESC/contigs.paths"

contig_graph = spades.get_contig_graph(graph_file, contig_paths_file)
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