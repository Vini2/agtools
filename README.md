# agtools: Tools for manipulating assembly graphs

This repo contains useful scripts to manipulate assembly graph files from different assemblers.

Currently, I have added code for [SPAdes](https://github.com/ablab/spades) (including metaSPAdes) assembly graphs.

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
git clone https://github.com/metagentools/gbintk.git

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

