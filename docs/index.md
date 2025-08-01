# agtools: Tools for manipulating assembly graphs

agtools is a toolkit for manipulating assembly graphs, with a focus on the [Graphical Fragment Assembly (GFA) format](https://github.com/GFA-spec/GFA-spec). It offers a command-line interface for tasks such as graph format conversion, segment filtering, and component extraction. Supported formats include [GFA](https://github.com/pmelsted/GFA-spec/blob/master/GFA-spec.md), [FASTG](https://web.archive.org/web/20211209213905/http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf), [ASQG](https://github.com/jts/sga/wiki/ASQG-Format) and [GraphViz DOT](http://www.graphviz.org/content/dot-language). Additionally, it provides a Python package interface that exposes assembler-specific functionality for advanced analysis and integration based on the GFA format.

## Quick install

Install using pip:

```shell
pip install agtools
```

Install using conda:

```shell
conda install -c bioconda agtools
```

Further details are available in the [Installation Guide](install.md).

## Documentation

**Tutorials**

* [CLI examples](examples.md)
* [API tutorial](tutorial.md)

### References

* [CLI reference](cli.md)
* [API reference](api.md)
* [Source code](https://github.com/Vini2/agtools)

### Versions and support

* [Changelog](changelog.md)
* [FAQ](faq.md)