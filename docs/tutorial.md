# API Tutorial

This page is a detailed tutorial of *agtools*' API. If you have not installed *agtools* yet, refer to [Installing *agtools*](install.md).

## Importing *agtools*

You can import *agtools* within a Python environment.

```python
% python
Python 3.13.5 | packaged by Anaconda, Inc. | (main, Jun 12 2025, 11:23:37) [Clang 14.0.6 ] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> import agtools
>>> print(agtools.__version__)
0.1.0
```

*agtools* provides two graph classes: `UnitigGraph` and `ContigGraph`. You can import them as follows.

```python
>>> from agtools.core.unitig_graph import UnitigGraph
>>> from agtools.core.contig_graph import ContigGraph
```

!!! note
    `UnitigGraph` and `ContigGraph` have slightly different implementations. If you want to just load a GFA file use the `UnitigGraph` class. Functions to load assembler-specific contig graphs are provided separately.

## Loading a GFA file

You can load a GFA file using the `from_gfa` method of the `UnitigGraph` class.

```python
ug = UnitigGraph.from_gfa("assembly_graph.gfa")
```

This will load the original segments (denoted by `S` tags) as vertices and links (denoted by `L` tags) as edges. If you are not familiar with the GFA format please refer to the [GFA Format Specification](https://gfa-spec.github.io/GFA-spec/GFA1.html).

You can view the different attributes of the graph.

```python
>>> ug.file_path
'assembly_graph.gfa'
>>> ug.vcount
34
>>> ug.ecount
23
```

You can call different functions to calculate graph and sequence based statistics.

```python
>>> ug.calculate_average_node_degree()
2
>>> ug.calculate_average_segment_length()
3494
>>> ug.calculate_n50_l50()
(15000, 12)
```

You can retrieve a sequence given the segment ID as follows.

```python
>>> ug.get_segment_sequence("unitig_1")
Seq('ATGCGTACGGGGTAAGTGAGCCTG')
```

!!! note
    Assembly graphs can be huge (10-100 GB in size). Hence, segment sequences are not loaded in to memory when creating the graph object. Instead, file pointers are kept for quick retrieval of sequences when needed. 


## Loading graphs from different assemblers

Different assemblers have different ways of representing assembly graphs. Some assemblers generated a unitig graph and resolved contigs from it where as some assemblers directly generate a contig graph. *agtools* currently supports two short-read assemblers SPAdes and MEGAHIT, and two long-read assemblers Flye and myloasm. 

### Loading a SPAdes graph

SPAdes generates a unitig graph and resolves longer paths as contigs. You will need three files to load a SPAdes contig graph.

* `assembly_graph_with_scaffolds.gfa`
* `contigs.fasta`
* `contigs.paths`

You can load a SPAdes contig graph as follows.

```python
>>> from agtools.assemblers import spades
>>> graph_file = "tests/data/ESC/assembly_graph_with_scaffolds.gfa"
>>> contigs_file = "tests/data/ESC/contigs.fasta"
>>> contig_paths_file = "tests/data/ESC/contigs.paths"
>>> contig_graph = spades.get_contig_graph(graph_file, contigs_file, contig_paths_file)
```

!!! note
    If you want to load the SPAdes unitig graph, you can load it as

    ```python
    unitig_graph = UnitigGraph.from_gfa("tests/data/ESC/assembly_graph_with_scaffolds.gfa")
    ```

### Loading a MEGAHIT graph

!!! note
    By default, MEGAHIT will not generate the assembly graph. You have to run `megahit_core` to get the assembly graph which will be in the FASTG format. Please refer to [MEGAHIT basic usage](https://github.com/voutcn/megahit?tab=readme-ov-file#basic-usage).

    ```bash
    megahit_core contig2fastg 141 out/intermediate_contigs/final.contigs.fa > final.fastg
    ```

    Then you can use the *agtools* CLI command `fastg2gfa` to convert it to GFA format. Please refer to the [CLI reference](cli.md) for further details on the `fastg2gfa` subcommand.

Once you have the GFA file of the assembly graph and the contigs file, you can load the contig graph as follows.

```python
>>> from agtools.assemblers import megahit
>>> graph_file = "tests/data/5G/final.gfa"
>>> contig_file = "tests/data/5G/final.contigs.fa"
>>> contig_graph = megahit.get_contig_graph(graph_file, contig_file)
```

### Loading a Flye graph

Similar to SPAdes, Flye also generates a unitig graph and resolves longer paths as contigs. You will need three files to load a Flye contig graph.

* `assembly_graph.gfa`
* `assembly.fasta`
* `assembly_info.txt`

You can load a Flye contig graph as follows.

```python
>>> from agtools.assemblers import flye
>>> graph_file = "tests/data/1Y3B/assembly_graph.gfa"
>>> contigs_file = "tests/data/1Y3B/assembly.fasta"
>>> contig_paths_file = "tests/data/1Y3B/assembly_info.txt"
>>> contig_graph = flye.get_contig_graph(graph_file, contigs_file, contig_paths_file)
```

!!! note
    If you want to load the Flye unitig graph, you can load it as

    ```python
    unitig_graph = UnitigGraph.from_gfa("tests/data/1Y3B/assembly_graph.gfa")
    ```

### Loading a myloasm graph

Myloasm generates an assembly graph in GFA format with contigs. However, not all the contigs that appear in the graph will appear in the FASTA file. 

!!! note
    Note that some contigs present in the myloasm assembly graph may not be included in the corresponding FASTA file. In such cases, you can use the *agtools* CLI command `clean` to GFA the file with contigs from the FASTA file. Please refer to the [CLI reference](cli.md) for further details on the `clean` subcommand.

You can load a myloasm contig graph as follows.

```python
>>> from agtools.assemblers import myloasm
>>> graph_file = "tests/data/myloasm/final_contig_graph.gfa"
>>> contigs_file = "tests/data/myloasm/assembly_primary.fa"
>>> contig_graph = myloasm.get_contig_graph(graph_file, contigs_file)
```