# **agtools**: A Software Framework to Manipulate Assembly Graphs

[![DOI](https://img.shields.io/badge/DOI-10.1093/bioadv/vbag126-blue)](https://doi.org/10.1093/bioadv/vbag126)
![GitHub License](https://img.shields.io/github/license/Vini2/agtools)
[![Python](https://img.shields.io/badge/Python-3776AB?logo=python&logoColor=fff)](#)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/agtools)
[![Conda](https://img.shields.io/conda/dn/bioconda/agtools)](https://anaconda.org/bioconda/agtools)
[![PyPI version](https://badge.fury.io/py/agtools.svg)](https://badge.fury.io/py/agtools)
[![Downloads](https://static.pepy.tech/badge/agtools)](https://pepy.tech/project/agtools)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1cO6gDsGqCIJUvnDcwpy-EuB3FoXG4M-R?usp=sharing)

[![CI](https://github.com/Vini2/agtools/actions/workflows/testing_python_app.yml/badge.svg)](https://github.com/Vini2/agtools/actions/workflows/testing_python_app.yml)
[![codecov](https://codecov.io/gh/Vini2/agtools/graph/badge.svg?token=nYzx0Pd0h6)](https://codecov.io/gh/Vini2/agtools)
[![CodeQL Advanced](https://github.com/Vini2/agtools/actions/workflows/codeql.yml/badge.svg)](https://github.com/Vini2/agtools/actions/workflows/codeql.yml)
![GitHub last commit (branch)](https://img.shields.io/github/last-commit/Vini2/agtools/main?color=8a35da)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Documentation Status](https://readthedocs.org/projects/agtools/badge/?version=latest)](https://agtools.readthedocs.io/en/latest/?badge=latest)

`agtools` (pronounced as *a-g-tools*) is a Python framework for manipulating assembly graphs for downstream metagenomic applications, with a focus on the [Graphical Fragment Assembly (GFA) format](https://github.com/GFA-spec/GFA-spec). It offers a command-line interface for tasks such as graph format conversion, segment filtering, and component extraction. Supported formats include [GFA](https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md), [FASTG](https://web.archive.org/web/20211209213905/http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf), [ASQG](https://github.com/jts/sga/wiki/ASQG-Format) and [GraphViz DOT](http://www.graphviz.org/content/dot-language). Additionally, it provides a Python package interface that exposes assembler-specific functionality for advanced analysis and integration based on the GFA format.

For detailed instructions on installation and usage, please refer to the [**documentation hosted at Read the Docs**](https://agtools.readthedocs.io).

## Installing `agtools`

### Using pip

You can install `agtools` from [PyPI](https://pypi.org/project/agtools/) using `pip`.

```bash
pip install agtools
```

### Using conda

You can install `agtools` from [Bioconda](https://anaconda.org/bioconda/agtools) using `conda` or [`mamba`](https://mamba.readthedocs.io/en/latest/index.html).

```bash
mamba install -c bioconda agtools
```

### For development

Please follow the steps below to install `agtools` using [`flit`](https://flit.pypa.io/en/stable/) for development.

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

  agtools: A Software Framework to Manipulate Assembly Graphs

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
  gfa2fastg  Convert GFA file to FASTG format
  asqg2gfa   Convert ASQG file to GFA format
  gfa2asqg   Convert GFA file to ASQG format
  gfa2dot    Convert GFA file to DOT format (GraphViz)
  gfa2fasta  Get segments in FASTA format
  gfa2adj    Get adjacency matrix of the assembly graph
```

## Documentation

Please refer to the complete documentation available at [Read the docs](https://agtools.readthedocs.io/).

## Issues and Questions

If you want to test (or break) `agtools` give it a try and report any issues and suggestions under [agtools Issues](https://github.com/Vini2/agtools/issues).

If you come across any questions, please have a look at the [agtools FAQ page](https://agtools.readthedocs.io/en/latest/faq/). If your question is not here, feel free to post it under [agtools Issues](https://github.com/Vini2/agtools/issues).

If you want to request a feature, please post it under [agtools Issues](https://github.com/Vini2/agtools/issues) as a feature request.

## Citation

`agtools` is publlished in [Bioinformatics Advances](https://doi.org/10.1093/bioadv/vbag126) with DOI: [10.1093/bioadv/vbag126](https://doi.org/10.1093/bioadv/vbag126). If you use `agtools` in your work, please cite as follows.

> Vijini Mallawaarachchi, George Bouras, Ryan R Wick, Susanna R Grigson, Bhavya Papudeshi, Robert A Edwards, agtools: A Software Framework to Manipulate Assembly Graphs, Bioinformatics Advances, Volume 6, Issue 1, 2026, vbag126, https://doi.org/10.1093/bioadv/vbag126

```bibtex
@article{10.1093/bioadv/vbag126,
    author = {Mallawaarachchi, Vijini and Bouras, George and Wick, Ryan R and Grigson, Susanna R and Papudeshi, Bhavya and Edwards, Robert A},
    title = {agtools: a software framework to manipulate assembly graphs},
    journal = {Bioinformatics Advances},
    volume = {6},
    number = {1},
    pages = {vbag126},
    year = {2026},
    month = {01},
    abstract = {Assembly graphs are a fundamental data structure used by genome and metagenome assemblers to represent sequences and their overlap information, facilitating the assembler in constructing longer genomic fragments. Apart from their core use in assemblers, assembly graphs have become increasingly important in a range of downstream applications such as metagenomic binning, plasmid detection, viral genome resolution, and haplotype phasing. However, there is a need for a comprehensive tool that allows programmatic access to manipulate assembly graphs (e.g. parse, convert, filter, and analyze) across different assembly graph formats.Here we present agtools, an open-source Python framework to manipulate assembly graphs produced by commonly used assemblers. agtools provides a command-line interface for tasks such as assembly graph format conversion, segment filtering, and component extraction. It also exposes a Python package interface to load, query, and analyze assembly graphs from popular genome and metagenome assemblers. This enables streamlined assembly-graph-based analyses that can be integrated into other bioinformatics software and workflows.The source code of agtools is hosted on GitHub at https://github.com/Vini2/agtools and the documentation is available at https://agtools.readthedocs.io/. agtools can also be installed from Bioconda (https://anaconda.org/bioconda/agtools) and PyPI (https://pypi.org/project/agtools/).},
    issn = {2635-0041},
    doi = {10.1093/bioadv/vbag126},
    url = {https://doi.org/10.1093/bioadv/vbag126},
    eprint = {https://academic.oup.com/bioinformaticsadvances/article-pdf/6/1/vbag126/68223878/vbag126.pdf},
}
```
