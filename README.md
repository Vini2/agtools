# agtools: Tools for manipulating assembly graphs

`agtools` is a toolkit for manipulating assembly graphs, with a focus on metagenomic applications. It offers a command-line interface for tasks such as graph format conversion, segment filtering, and subgraph extraction. Additionally, it provides a Python package interface that exposes assembler-specific functionality for advanced analysis and integration.

## Requirements

You should have Python and the following packages installed.

* [flit](https://flit.pypa.io/en/stable/)
* [click](https://click.palletsprojects.com/en/stable/)
* [python-igraph](https://python.igraph.org/en/stable/index.html)
* [biopython](https://biopython.org/)
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
  filter     Filter segments from GFA file
  component  Extract a component containing a given segment
  fastg2gfa  Convert FASTG file to GFA format
  gfa2fastg  Convert GFA file to FASTG format
  gfa2dot    Convert GFA file to DOT format (Graphviz)
  gfa2fasta  Get segments in FASTA format
```