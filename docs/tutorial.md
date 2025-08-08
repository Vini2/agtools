# API Tutorial

This page is a detailed tutorial of *agtools*' API. If you want to get a quick idea on how the *agtools* API works, do check out the [Quick Start Guide](quickstart.md) If you have not installed *agtools* yet, refer to [Installing *agtools*](install.md).

## Importing *agtools*

You can import *agtools* within a Python environment.

```python
% python
Python 3.13.5 | packaged by Anaconda, Inc. | (main, Jun 12 2025, 11:23:37) [Clang 14.0.6 ] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> import agtools
>>> print(agtools.__version__)
0.1.2
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

You can retrieve a sequence given the segment ID as follows. A `Bio.Seq.Seq` object will be returned.

```python
>>> seq = ug.get_segment_sequence("unitig_1")
>>> seq
Seq('ATGCGTACGGGGTAAGTGAGCCTG')
>>> seq.reverse_complement()
Seq('CAGGCTCACTTACCCCGTACGCAT')
```

!!! note
    Assembly graphs can be huge (10-100 GB in size). Hence, segment sequences are not loaded in to memory when creating the graph object. Instead, file pointers are kept for quick retrieval of sequences when needed. 


## Loading contig graphs

Different assemblers have different ways of representing assembly graphs. Some assemblers generated a unitig graph and resolved contigs from it where as some assemblers directly generate a contig graph. *agtools* currently supports the following short-read and long-read assemblers. 

* [SPAdes](https://github.com/ablab/spades)
* [MEGAHIT](https://github.com/voutcn/megahit)
* [Flye](https://github.com/mikolmogorov/Flye)
* [myloasm](https://github.com/bluenote-1577/myloasm)

Please refer to the [assembler-specific examples](assemblerexamples.md) for further details on the graph representations and more details examples.

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
>>> cg = spades.get_contig_graph(graph_file, contigs_file, contig_paths_file)
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
>>> cg = megahit.get_contig_graph(graph_file, contig_file)
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
>>> cg = flye.get_contig_graph(graph_file, contigs_file, contig_paths_file)
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
>>> cg = myloasm.get_contig_graph(graph_file, contigs_file)
```

## Querying contig graphs

You can view the different attributes of the contig graphs obtained using the assembler-specific modules as follows.

```python
>>> cg.file_path
'tests/data/ESC/assembly_graph_with_scaffolds.gfa'
>>> cg.vcount
189
>>> cg.ecount
394
```

You can get the mapping of the internal node ID to the contig names as follows.

```python
>>> cg.contig_names
bidict({0: 'NODE_1_length_488682_cov_86.190505', 1: 'NODE_2_length_472233_cov_17.669606', 2: 'NODE_3_length_354360_cov_17.661738', 3: 'NODE_4_length_346431_cov_86.228266', ...})
```

!!! note
    The internal ID starting from 0 is used to index the nodes in the `igraph` object. The corresponding contig name can be obtained from the `contig_names` attribute of the contig graph object. This is useful when traversing nodes using the `igraph` object.

You can call different functions to calculate graph and sequence based statistics.

```python
>>> cg.calculate_average_node_degree()
4
>>> cg.calculate_total_length()
8341464
>>> cg.calculate_average_contig_length()
44134
>>> cg.calculate_n50_l50()
(220639, 14)
>>> cg.get_gc_content()
0.4350709899365387
```

You can retrieve a sequence given the segment ID as follows. A Bio.Seq.Seq object will be returned.

```python
>>> cg.get_contig_sequence("NODE_189_length_56_cov_33.000000")
Seq('TGGCTCTTCAGGATCCAGGGTGTAGTCGGGGTCTGAATCCTCCGGTCTCCAGGAGG')
```

You can get the neighbours of a contig as follows.

```python
>>> cg.get_neighbors("NODE_1_length_488682_cov_86.190505")
['NODE_4_length_346431_cov_86.228266', 'NODE_9_length_265823_cov_88.260370', 'NODE_14_length_220639_cov_88.091915', 'NODE_19_length_193605_cov_88.758791', 'NODE_21_length_158927_cov_87.573997', 'NODE_23_length_117248_cov_87.447902', 'NODE_29_length_97556_cov_86.828022', 'NODE_33_length_87414_cov_86.098490', 'NODE_42_length_49567_cov_91.112902', 'NODE_44_length_45842_cov_86.030074', 'NODE_46_length_42994_cov_88.472018', 'NODE_47_length_42793_cov_92.771094', 'NODE_50_length_41992_cov_86.797863', 'NODE_63_length_17530_cov_92.904492', 'NODE_89_length_395_cov_158.788235', 'NODE_95_length_284_cov_84.877729', 'NODE_108_length_160_cov_106.523810', 'NODE_118_length_129_cov_90.905405', 'NODE_136_length_114_cov_87.813559', 'NODE_139_length_111_cov_73.875000', 'NODE_146_length_99_cov_86.818182', 'NODE_148_length_98_cov_74.906977', 'NODE_154_length_88_cov_148.030303', 'NODE_156_length_85_cov_89.700000', 'NODE_164_length_65_cov_81.100000', 'NODE_167_length_63_cov_149.000000']
```

You can check if two contigs are connected by a path as follows.

```python
>>> cg.is_connected("NODE_1_length_488682_cov_86.190505", "NODE_146_length_99_cov_86.818182")
True
```