# API Tutorial

This page is a detailed tutorial of *agtools*' API. If you want to get a quick overview of how the *agtools* API works, please refer to the [Quick Start Guide](quickstart.md). If you have not installed *agtools* yet, see [Installing *agtools*](install.md).

## Importing *agtools*

You can import *agtools* within a Python environment.

```python
% python
Python 3.13.5 | packaged by Anaconda, Inc. | (main, Jun 12 2025, 11:23:37) [Clang 14.0.6 ] on darwin
Type "help", "copyright", "credits" or "license" for more information.
>>> import agtools
>>> print(agtools.__version__)
1.0.2
```

*agtools* provides two graph classes: `UnitigGraph` and `ContigGraph`. You can import them as follows.

```python
>>> from agtools.core.unitig_graph import UnitigGraph
>>> from agtools.core.contig_graph import ContigGraph
```

!!! note
    `UnitigGraph` and `ContigGraph` have slightly different implementations. If you only want to load a GFA file, use the `UnitigGraph` class. Functions to load assembler-specific contig graphs are provided separately.

## Loading a GFA file

You can load a GFA file using the `from_gfa` method of the `UnitigGraph` class. Let's load an example [GFA file](https://github.com/Vini2/agtools/tree/main/tests/data/ESC).

```python
ug = UnitigGraph.from_gfa("assembly_graph_with_scaffolds.gfa")
```

This will load the original segments (denoted by `S` tags) as vertices and links (denoted by `L` tags) as edges. If you are not familiar with the GFA format, please refer to the [GFA Format Specification](https://gfa-spec.github.io/GFA-spec/GFA1.html).

You can view the different attributes of the graph.

```python
>>> ug.file_path
'assembly_graph_with_scaffolds.gfa'
>>> ug.vcount   # number of vertices (segments) in the graph
982
>>> ug.ecount   # number of edges in the graph
1265
>>> ug.lcount   # number of lines starting with tag 'L'
1318
>>> ug.pcount   # number of lines starting with tag 'P'
190
```

You can call various functions to calculate graph- and sequence-based statistics.

```python
>>> ug.calculate_average_node_degree()
2.576374745417515
>>> ug.calculate_average_segment_length()
8490.319755600814
>>> ug.calculate_n50_l50()
(60706, 36)
```

You can retrieve a sequence given the segment ID as follows. A `Bio.Seq.Seq` object will be returned.

```python
>>> seq = ug.get_segment_sequence("3042")
>>> seq
Seq('TGATTTTCGCGCGATTACTACGATGATTTCAAACGATTCCTCTGATTATTTCACGC')
>>> seq.reverse_complement()
Seq('GCGTGAAATAATCAGAGGAATCGTTTGAAATCATCGTAGTAATCGCGCGAAAATCA')
```

!!! note
    Assembly graphs can be very large (10–100 GB). Therefore, segment sequences are not loaded into memory when creating the graph object. Instead, file pointers are stored for efficient sequence retrieval.


## Loading contig graphs

Different assemblers represent assembly graphs in different ways. Some assemblers generate a unitig graph and resolve contigs from it, whereas others directly generate contig graphs. *agtools* currently supports the following assemblers:

* [SPAdes](https://github.com/ablab/spades)
* [MEGAHIT](https://github.com/voutcn/megahit)
* [Flye](https://github.com/mikolmogorov/Flye)
* [myloasm](https://github.com/bluenote-1577/myloasm)

Please refer to the [assembler-specific examples](assemblerexamples.md) for further details on the graph representations and more detailed examples.

### Loading a SPAdes graph

SPAdes generates a unitig graph and resolves longer paths as contigs. You will need three files to load a SPAdes contig graph:

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
    If you only want to load the SPAdes unitig graph, you can do so using:

    ```python
    unitig_graph = UnitigGraph.from_gfa("tests/data/ESC/assembly_graph_with_scaffolds.gfa")
    ```

### Loading a MEGAHIT graph

!!! note
    By default, MEGAHIT does not generate the assembly graph. You must run `megahit_core` to produce it in FASTG format. See [MEGAHIT basic usage](https://github.com/voutcn/megahit?tab=readme-ov-file#basic-usage).

    ```bash
    megahit_core contig2fastg 141 out/intermediate_contigs/final.contigs.fa > final.fastg
    ```

    You can then convert it to GFA format using the *agtools* CLI command `fastg2gfa`. See the [CLI reference](cli.md) for details.

Once you have the GFA file and the contigs file, you can load the contig graph as follows.

```python
>>> from agtools.assemblers import megahit
>>> graph_file = "tests/data/5G/final.gfa"
>>> contig_file = "tests/data/5G/final.contigs.fa"
>>> cg = megahit.get_contig_graph(graph_file, contig_file)
```

### Loading a Flye graph

Similar to SPAdes, Flye generates a unitig graph and resolves contigs from it. You will need the following files:

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

!!! warning
    The final Flye contig sequences may differ from those derived from the GFA file due to polishing steps and the way Flye constructs contigs (see [Flye/issues/610](https://github.com/mikolmogorov/Flye/issues/610)).

!!! note
    To load only the Flye unitig graph:

    ```python
    unitig_graph = UnitigGraph.from_gfa("tests/data/1Y3B/assembly_graph.gfa")
    ```


### Loading a myloasm graph

Myloasm generates a GFA assembly graph with contigs. However, not all contigs present in the graph necessarily appear in the FASTA file.

!!! note
    Some contigs in the myloasm GFA may be missing from the FASTA file. In such cases, you can use the *agtools* CLI command `clean` to filter the GFA file using the FASTA file. See the [CLI reference](cli.md) for details.

You can load a myloasm contig graph as follows.

```python
>>> from agtools.assemblers import myloasm
>>> graph_file = "tests/data/myloasm/final_contig_graph.gfa"
>>> contigs_file = "tests/data/myloasm/assembly_primary.fa"
>>> cg = myloasm.get_contig_graph(graph_file, contigs_file)
```

## Querying contig graphs

You can inspect the properties of a contig graph obtained using the assembler-specific modules as follows.

```python
>>> cg.file_path
'tests/data/ESC/assembly_graph_with_scaffolds.gfa'
>>> cg.vcount
189
>>> cg.ecount
394
>>> cg.lcount
1318
```

You can retrieve contig names as follows.

```python
>>> cg.contig_names
['NODE_1_length_488682_cov_86.190505', 'NODE_2_length_472233_cov_17.669606', 'NODE_3_length_354360_cov_17.661738', ...]
```

You can obtain the mapping from contig name to internal ID:

```python
>>> cg.contig_name_to_id
{'NODE_1_length_488682_cov_86.190505': 0, 'NODE_2_length_472233_cov_17.669606': 1, 'NODE_3_length_354360_cov_17.661738': 2, ...}
```

!!! note
    Internal IDs start from 0 and are used to index nodes in the underlying `igraph` object. The internal ID of each contig can be obtained from `contig_name_to_id`. The corresponding contig name can be obtained from `contig_names`. These mappings are useful when traversing nodes using the `igraph` object.

You can calculate graph and sequence statistics as follows.

```python
>>> cg.calculate_average_node_degree()
4.169312169312169
>>> cg.calculate_total_length()
8341464
>>> cg.calculate_average_contig_length()
44134.730158730155
>>> cg.calculate_n50_l50()
(220639, 14)
>>> cg.get_gc_content()
0.4350709899365387
```

You can retrieve a sequence given the segment ID as follows. A `Bio.Seq.Seq` object will be returned.

```python
>>> cg.get_contig_sequence("NODE_189_length_56_cov_33.000000")
Seq('TGGCTCTTCAGGATCCAGGGTGTAGTCGGGGTCTGAATCCTCCGGTCTCCAGGAGG')
```

You can retrieve neighbouring contigs as follows.

```python
>>> cg.get_neighbors("NODE_1_length_488682_cov_86.190505")
['NODE_4_length_346431_cov_86.228266', 'NODE_9_length_265823_cov_88.260370', 'NODE_14_length_220639_cov_88.091915', 'NODE_19_length_193605_cov_88.758791', 'NODE_21_length_158927_cov_87.573997', 'NODE_23_length_117248_cov_87.447902', 'NODE_29_length_97556_cov_86.828022', 'NODE_33_length_87414_cov_86.098490', 'NODE_42_length_49567_cov_91.112902', 'NODE_44_length_45842_cov_86.030074', 'NODE_46_length_42994_cov_88.472018', 'NODE_47_length_42793_cov_92.771094', 'NODE_50_length_41992_cov_86.797863', 'NODE_63_length_17530_cov_92.904492', 'NODE_89_length_395_cov_158.788235', 'NODE_95_length_284_cov_84.877729', 'NODE_108_length_160_cov_106.523810', 'NODE_118_length_129_cov_90.905405', 'NODE_136_length_114_cov_87.813559', 'NODE_139_length_111_cov_73.875000', 'NODE_146_length_99_cov_86.818182', 'NODE_148_length_98_cov_74.906977', 'NODE_154_length_88_cov_148.030303', 'NODE_156_length_85_cov_89.700000', 'NODE_164_length_65_cov_81.100000', 'NODE_167_length_63_cov_149.000000']
```

You can check whether two contigs are connected as follows.

```python
>>> cg.is_connected("NODE_1_length_488682_cov_86.190505", "NODE_146_length_99_cov_86.818182")
True
```