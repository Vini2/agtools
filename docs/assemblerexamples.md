# Assembler-specific examples

This page guides you through the assembler-specific contig graph representations provided by the *agtools* API, along with explanations of how contig graphs are constructed.

---

## SPAdes contig graph

Files required:

* GFA file - `assembly_graph_with_scaffolds.gfa`
* Contigs file - `contigs.fasta`
* Contig paths file - `contigs.paths`

Below is an example SPAdes assembly graph (note that this is not a real SPAdes output).

```text
H	VN:Z:1.0
S	8925540	ACGTAC
S	8925544	CGTACGT
S	8925548	GTACGTA
S	8925552	TACGTAG
L	8925540	+	8925544	+	5M
L	8925544	+	8925548	+	5M
L	8925548	+	8925552	+	5M
P	NODE_1_length_13_cov_20.0	8925540+,8925544+	5M,5M
P	NODE_2_length_21_cov_18.5	8925544+,8925548+,8925552+	5M,5M,5M
```

SPAdes generates a unitig graph and resolves longer paths as contigs. *agtools* constructs the contig graph based on prefix–suffix matching of contigs. In other words, if the last segment of one contig matches the first segment of another, a link is created.

For example, the path `NODE_1_length_13_cov_20.0` ends with segment `8925544`, and `NODE_2_length_21_cov_18.5` starts with `8925544`.


You can load a SPAdes contig graph as follows. We will use the provided [test data](https://github.com/Vini2/agtools/tree/main/tests/data/ESC).

```python
>>> from agtools.assemblers import spades
>>> graph_file = "tests/data/ESC/assembly_graph_with_scaffolds.gfa"
>>> contigs_file = "tests/data/ESC/contigs.fasta"
>>> contig_paths_file = "tests/data/ESC/contigs.paths"
>>> cg = spades.get_contig_graph(graph_file, contigs_file, contig_paths_file)
```

You can check the number of vertices (contigs) and edges (connections):

```python
>>> cg.vcount
189
>>> cg.ecount
394
```

!!! note
    There are fewer contigs than segments (`S` tags) in the graph (189 contigs and 982 segments in this example) because SPAdes resolves paths with multiple segments into longer contigs. Each contig is made up of one or more segments.

!!! note
    There are also fewer edges than links (`L` tags) in the graph (394 edges and 1318 links in this example). This is because *agtools* only represents prefix-suffix overlaps of contigs. Also, the graph is undirected and simplified to remove multiple edges and self-loops.


## MEGAHIT contig graph

MEGAHIT assembly graphs represent contig sequences as segments and overlaps between contigs as links, so they are straightforward to interpret. However, when building the graph file from `megahit_core`, the contigs are renamed.

For example, the first few entries of `final.contigs.fa` look like:

```text
>k141_4704 flag=0 multi=1.0000 len=301
TTCTGCAACACTTAAAGACATTCAAACTTTTTATGTTTGGGGTTACTGTAAAAACGATCAAATTAAGTGGTCCGTAGCCAAGGGTTTAATCGATAAAGAAGAATACACTTTGATCACTGGAGAAAAATATCCAGAAACAAAAGATGAAAAGTCACAGGTGTAATACTTGTGGCTTTTTAATTTGAATAAAGTGGGTGGCATAATGTTTGGATTTACCAAACGACATGAACAAGATTGGAGTTTAACGCGATTAGAAGAAAATGATAAGACTATGTTTGAAAAATTCGACAGAATAGAAGAT
>k141_10879 flag=1 multi=1.0000 len=301
GCAGTTAAGGCAGATCAAAATTTCAAGCGTACAGATGATAAGACTAAAGATGTTGAACTAAGTGATAATGATAAAAAATTACAAAAAACTATCAAGGTTGCTAAGAAGACTTTAGATTCGTTCGATAATGAAAAAGCAAAAGCTAAATTAGATGCTAAAATAGAAGATTTGAAACAAAAAGTATTAGAAGCAAGTTTTGAATTAAATCAATAAGATTCAAAAGAAGTTACACCAGAAGTTAAGTTAGAAAAACAAAAGTTAATTAAAGATATCACTGAAAAAGAAGCTAAGTTATCCGAAC
...
```

The `final.graph.fastg` file contains entries such as:

```text
>NODE_1_length_205_cov_1.0000_ID_1;
TTCTGCAACACTTAAAGACATTCAAACTTTTTATGTTTGGGGTTACTGTAAAAACGATCAAATTAAGTGGTCCGTAGCCAAGGGTTTAATCGATAAAGAAGAATACACTTTGATCACTGGAGAAAAATATCCAGAAACAAAAGATGAAAAGTCACAGGTGTAATACTTGTGGCTTTTTAATTTGAATAAAGTGGGTGGCATAATGTTTGGATTTACCAAACGACATGAACAAGATTGGAGTTTAACGCGATTAGAAGAAAATGATAAGACTATGTTTGAAAAATTCGACAGAATAGAAGAT
>NODE_1_length_205_cov_1.0000_ID_1';
ATCTTCTATTCTGTCGAATTTTTCAAACATAGTCTTATCATTTTCTTCTAATCGCGTTAAACTCCAATCTTGTTCATGTCGTTTGGTAAATCCAAACATTATGCCACCCACTTTATTCAAATTAAAAAGCCACAAGTATTACACCTGTGACTTTTCATCTTTTGTTTCTGGATATTTTTCTCCAGTGATCAAAGTGTATTCTTCTTTATCGATTAAACCCTTGGCTACGGACCACTTAATTTGATCGTTTTTACAGTAACCCCAAACATAAAAAGTTTGAATGTCTTTAAGTGTTGCAGAA
>NODE_2_length_301_cov_1.0000_ID_3;
GCAGTTAAGGCAGATCAAAATTTCAAGCGTACAGATGATAAGACTAAAGATGTTGAACTAAGTGATAATGATAAAAAATTACAAAAAACTATCAAGGTTGCTAAGAAGACTTTAGATTCGTTCGATAATGAAAAAGCAAAAGCTAAATTAGATGCTAAAATAGAAGATTTGAAACAAAAAGTATTAGAAGCAAGTTTTGAATTAAATCAATAAGATTCAAAAGAAGTTACACCAGAAGTTAAGTTAGAAAAACAAAAGTTAATTAAAGATATCACTGAAAAAGAAGCTAAGTTATCCGAAC
>NODE_2_length_301_cov_1.0000_ID_3';
GTTCGGATAACTTAGCTTCTTTTTCAGTGATATCTTTAATTAACTTTTGTTTTTCTAACTTAACTTCTGGTGTAACTTCTTTTGAATCTTATTGATTTAATTCAAAACTTGCTTCTAATACTTTTTGTTTCAAATCTTCTATTTTAGCATCTAATTTAGCTTTTGCTTTTTCATTATCGAACGAATCTAAAGTCTTCTTAGCAACCTTGATAGTTTTTTGTAATTTTTTATCATTATCACTTAGTTCAACATCTTTAGTCTTATCATCTGTACGCTTGAAATTTTGATCTGCCTTAACTGC
...
```

Here

* `k141_4704` corresponds to `NODE_1_length_205_cov_1.0000_ID_1`
* `k141_10879` corresponds to `NODE_2_length_301_cov_1.0000_ID_3`

The `megahit` module of *agtools* automatically resolves this mapping.

Files required:

* GFA file - You can convert the `final.graph.fastg` in to GFA format using the *agtools* [`fastg2gfa`](https://agtools.readthedocs.io/en/latest/cli/#fastg2gfa) subcommand.
* Contigs file - `final.contigs.fa`

Once you convert the FASTG file to GFA format, you can load a MEGAHIT contig graph as follows. We will use the provided [test data](https://github.com/Vini2/agtools/tree/main/tests/data/ESC).

```python
>>> from agtools.assemblers import megahit
>>> graph_file = "tests/data/5G/final.gfa"
>>> contig_file = "tests/data/5G/final.contigs.fa"
>>> cg = megahit.get_contig_graph(graph_file, contig_file)
```

You can get the contig name mapping between the graph file and the FASTA file as follows.

```python
>>> cg.graph_to_contig_map
bidict({'NODE_1_length_205_cov_1.0000_ID_1': 'k141_4704', 'NODE_2_length_301_cov_1.0000_ID_3': 'k141_10879', 'NODE_3_length_301_cov_1.0000_ID_5': 'k141_9891', ...})
```

Given the contig name in the graph, you can get the contig name in the FASTA file as follows.

```python
>>> cg.graph_to_contig_map['NODE_1_length_301_cov_1.0000_ID_1']
'k141_4704'
```

You can query the other way around as well. Given the contig name in the FASTA file, you can get the contig name in the graph as follows.

```python
>>> cg.graph_to_contig_map.inverse['k141_4704']
['NODE_1_length_205_cov_1.0000_ID_1']
```

!!! note 
    `graph_to_contig_map` is a bidirectional dictionary. Hence, you can query both ways: contig names from graph to FASTA and from FASTA to graph.

You can get the contig sequences using the contig name in the graph as follows.

```python
>>> cg.get_contig_sequence('NODE_1_length_205_cov_1.0000_ID_1')
Seq('CGCCCTCGCTCATGCACGATTGCGAAGGCTCCGGCATGCCATCTGAAGAAACAG...GAA')
```

Note that the contigs have some descriptions after the name in the FASTA file. You can get them as follows.

```python
>>> cg.contig_descriptions['k141_4704']
'k141_4704 flag=0 multi=1.0000 len=205'
```

## Flye contig graph

Flye contig graphs follow a similar structure to SPAdes, but links represent direct continuations rather than overlaps. Similar to SPAdes, *agtools* builds the contig graph based on prefix and suffix matching of contigs, *i.e.*, if the last segment of the first contig matches the first segment of the second contig, then a link between the two contigs is created.

Files required:

* GFA file - `assembly_graph.gfa`
* Contigs file - `assembly.fasta`
* Contig paths file - `assembly_info.txt`

You can load a Flye contig graph as follows. We will use the provided [test data](https://github.com/Vini2/agtools/tree/main/tests/data/1Y3B).

```python
>>> from agtools.assemblers import flye
>>> graph_file = "tests/data/1Y3B/assembly_graph.gfa"
>>> contigs_file = "tests/data/1Y3B/assembly.fasta"
>>> contig_paths_file = "tests/data/1Y3B/assembly_info.txt"
>>> cg = flye.get_contig_graph(graph_file, contigs_file, contig_paths_file)
```

You can check the number of vertices (contigs), edges (connections) and number of link lines in the GFA file.

```python
>>> cg.vcount
67
>>> cg.ecount
2
>>> cg.lcount
10
```

!!! note
    There are fewer contigs than segments (`S` tags) in the graph (67 contigs and 69 segments in this example). This is because Flye resolves longer segment paths as contigs. Each contig is made up of one or more segments.

!!! note
    There are fewer edges than links (`L` tags) in the graph (2 edges and 10 links in this example). This is because *agtools* only represents prefix-suffix overlaps of contigs. Also, the graph is undirected and simplified to remove duplicate edges.


## myloasm contig graph

myloasm outputs contig graphs directly, where segments represent contigs and links represent overlaps, so they are straightforward to represent.

Files required:

* GFA file - `final_contig_graph.gfa`
* Contigs file - `assembly_primary.fa`


```python
>>> from agtools.assemblers import myloasm
>>> graph_file = "tests/data/myloasm/final_contig_graph.gfa"
>>> contigs_file = "tests/data/myloasm/assembly_primary.fa"
>>> cg = myloasm.get_contig_graph(graph_file, contigs_file)
```

!!! warning
    Not all contigs in the `final_contig_graph.gfa` file appear in `assembly_primary.fa`, because contigs that are similar to larger contigs are removed (see [GitHub issue](https://github.com/bluenote-1577/myloasm/issues/5)). Also, some contig sequences are not present in the GFA file. This can raise warnings when retrieving contig sequences. Hence, it is recommended to run `agtools clean` to remove missing contigs and add back contig sequences to the GFA file.


You can check the number of vertices (contigs), edges (connections), the number of link lines and contig names:

```python
>>> cg.vcount
2
>>> cg.ecount
1
>>> cg.lcount
2
>>> cg.contig_names
['u913838ctg', 'u579439ctg']
>>> cg.contig_name_to_id
{'u913838ctg': 0, 'u579439ctg': 1}
```
