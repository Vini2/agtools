# agtools: A Software Framework to Manipulate Assembly Graphs

*agtools* is a Python framework for manipulating assembly graphs for downstream metagenomic applications, with a focus on the [Graphical Fragment Assembly (GFA) format](https://github.com/GFA-spec/GFA-spec). It offers a command-line interface for tasks such as graph format conversion, segment filtering, and component extraction. Supported formats include [GFA](https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md), [FASTG](https://web.archive.org/web/20211209213905/http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf), [ASQG](https://github.com/jts/sga/wiki/ASQG-Format) and [GraphViz DOT](http://www.graphviz.org/content/dot-language). In addition, it provides a Python package interface that exposes assembler-specific functionality for advanced analysis and integration based on the GFA format.

## Quick install

Install using [pip](https://pypi.org/project/agtools/):

```shell
pip install agtools
```

Install from the [Bioconda](https://anaconda.org/bioconda/agtools) distribution using [`conda`](https://docs.conda.io/projects/conda/en/latest/index.html) or [`mamba`](https://mamba.readthedocs.io/en/latest/index.html):

```shell
mamba install -c bioconda agtools
```

Further details are available in the [Installation Guide](install.md).

## Documentation

**Tutorials**

* [CLI examples](examples/stats.md)
* [API tutorial](tutorial.md)
* [Assembler-specific examples](assemblerexamples.md)
* [More detailed examples](moreexamples.md)
* [Example applications](exampleapplications.md)

**References**

* [CLI reference](cli.md)
* [API reference](api.md)
* [File formats](fileformats.md)
* [Source code](https://github.com/Vini2/agtools)

**Support**

* [Changelog](changelog.md)
* [FAQ](faq.md)

## Citation

*agtools* is currently under review. In the meantime, if you use *agtools* in your work, please cite the preprint as follows.

> Vijini Mallawaarachchi, George Bouras, Ryan R. Wick, Susanna R. Grigson, Bhavya Papudeshi, Robert A. Edwards; agtools: a software framework to manipulate assembly graphs; bioRxiv 2025.09.14.676178; doi: https://doi.org/10.1101/2025.09.14.676178