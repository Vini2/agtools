# Installing *agtools*

It is recommended to install *agtools* using either [PyPI](https://pypi.org/project/agtools/) or [Conda](https://anaconda.org/bioconda/agtools).

## PyPI

To install *agtools* globally, use the following command.

```bash
pip install agtools
```

## Conda

You can install *agtools* through [Bioconda](https://anaconda.org/bioconda/agtools) using [`conda`](https://docs.conda.io/projects/conda/en/latest/index.html) or [`mamba`](https://mamba.readthedocs.io/en/latest/index.html).

```bash
mamba install -c bioconda agtools
```

If you prefer to install *agtools* in a dedicated environment, you can do so as follows.

```bash
mamba create -n agtools
mamba activate agtools
mamba install -c bioconda agtools
```

Alternatively:

```bash
mamba create -n agtools -c bioconda agtools
mamba activate agtools
```

## Installing *agtools* from source

If you want to use the development version of *agtools*, you can install it using `flit` as follows. Please ensure that [`flit`](https://flit.pypa.io/en/stable/) is installed beforehand.

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

If you prefer not to create a dedicated environment and only want to install *agtools*, ensure the following dependencies are installed:

* [`flit`](https://flit.pypa.io/en/stable/) - for installation
* [`click`](https://click.palletsprojects.com/en/stable/) - for CLI argument parsing
* [`loguru`](https://loguru.readthedocs.io/en/stable/) - for logging
* [`bidict`](https://bidict.readthedocs.io/en/main/intro.html) - for bidirectional lookup
* [`python-igraph`](https://python.igraph.org/en/stable/index.html) - for graph operations
* [`biopython`](https://biopython.org/) - for sequence operations
* [`pandas`](https://pandas.pydata.org/) - for dataframes
