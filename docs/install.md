# Installing *agtools*

It is recommended to install *agtools* using either [PyPI](https://pypi.org/project/agtools/) or Conda.

## PyPI

To install *agtools* globally, use the following command.

```bash
pip install agtools
```

## Conda

You can install *agtools* through Bioconda.

```bash
conda install -c bioconda agtools
```

If you prefer to install *agtools* in your own environment, you can do so as follows.

```bash
conda create -n agtools
conda activate agtools
conda install -c bioconda agtools
```

OR

```bash
conda create -n agtools -c bioconda agtools
conda activate agtools
```

## Installing *agtools* from source

If you want to use the development version of *agtools*, you can install it using `flit` as follows. Please make sure you have [`flit`](https://flit.pypa.io/en/stable/) installed beforehand.

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

If you don't want to create your own environment but just install *agtools*, please make sure you have the following packages installed.

* [flit](https://flit.pypa.io/en/stable/)
* [click](https://click.palletsprojects.com/en/stable/)
* [loguru](https://loguru.readthedocs.io/en/stable/)
* [bidict](https://bidict.readthedocs.io/en/main/intro.html)
* [python-igraph](https://python.igraph.org/en/stable/index.html)
* [biopython](https://biopython.org/)
* [pandas](https://pandas.pydata.org/)
