# Frequently Asked Questions

Below are answers to common questions about using `agtools` for working with assembly graphs.

## What is `agtools`?

`agtools` provides a command line interface to manipulate assembly graphs and a Python package interface to to explore assembly graphs in your own software. You can refer the following resources to get an idea of the capabilities of `agtools`.

* [CLI examples](examples/stats.md)
* [CLI reference](cli.md)
* [API tutorial](tutorial.md)
* [API reference](api.md)

## Why would I need `agtools`?

Assembly graphs are gaining popularity, especially in the metagenomics space as it holds valuable connectivity information of the assembled sequences. This information can be extremely useful in downstream applications such as 

* Metagenomic binning - GraphBin, MetaCoAG, GraphMB
* Resolving bacterial strains - STRONG
* Resolving viral strains - VStrains, Phables
* Identifying plasmids - GraphPlas
* Sequence classification - 3CAC, 4CAC

where many solutions that use assembly graphs have been published.

Many applications require a programmatic and modular way to explore and manipulate assembly graphs. However, no standardiszed solution currently exists and developers often resort to writing custom, one-off code. `agtools` bridges this gap by offering a reusable and extensible toolkit for working with assembly graphs in a consistent and accessible manner.