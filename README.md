# agtools: Tools for manipulating assembly graphs

![GitHub License](https://img.shields.io/github/license/Vini2/agtools)
[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](#)
[![CI](https://github.com/Vini2/agtools/actions/workflows/testing_python_app.yml/badge.svg)](https://github.com/Vini2/agtools/actions/workflows/testing_python_app.yml)
[![codecov](https://codecov.io/gh/Vini2/agtools/graph/badge.svg?token=nYzx0Pd0h6)](https://codecov.io/gh/Vini2/agtools)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/Vini2/agtools/main?color=8a35da)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

`agtools` is a toolkit for manipulating assembly graphs, with a focus on the [Graphical Fragment Assembly (GFA) format](https://github.com/GFA-spec/GFA-spec). It offers a command-line interface for tasks such as graph format conversion, segment filtering, and component extraction. Supported formats include [GFA](https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md), [FASTG](https://web.archive.org/web/20211209213905/http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf), [ASQG](https://github.com/jts/sga/wiki/ASQG-Format) and [GraphViz DOT](http://www.graphviz.org/content/dot-language). Additionally, it provides a Python package interface that exposes assembler-specific functionality for advanced analysis and integration based on the GFA format.

## Requirements

You should have Python and the following packages installed.

* [flit](https://flit.pypa.io/en/stable/)
* [click](https://click.palletsprojects.com/en/stable/)
* [loguru](https://loguru.readthedocs.io/en/stable/)
* [bidict](https://bidict.readthedocs.io/en/main/intro.html)
* [python-igraph](https://python.igraph.org/en/stable/index.html)
* [biopython](https://biopython.org/)
* [pandas](https://pandas.pydata.org/)

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

## Documentation

Please refer to the complete documentation available at [Read the docs](https://agtools.readthedocs.io/)