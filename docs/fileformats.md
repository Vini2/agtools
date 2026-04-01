# Assembly Graph File Formats

Assembly graphs describe the structure of DNA assemblies using nodes (segments or contigs) and edges (connections or overlaps). Several file formats are used to represent these graphs, each with its own syntax and conventions.

This page summarises the four most common formats supported by *agtools*.

## GFA

**Extension**: `.gfa`

**Versions supported**: GFA1 (commonly used in assembly)

GFA (Graphical Fragment Assembly) is a standard format for describing genome assemblies as graphs, where:

* Headers (start with the tag `H`) denote metadata.
* Segments (start with the tag `S`) are nodes (contigs or unitigs).
* Links (start with the tag `L`) define edges between segments.
* Paths (start with the tag `P`) optionally define traversal orders.

*agtools* also supports the tags `J`, `C`, and `W`. Please refer to the [GFA specification](https://gfa-spec.github.io/GFA-spec/GFA1.html) for further details.

**Example**

```text
S	unitig1	ATGCGTACGGGGTAAGTGAGCCTG
S	unitig2	CGGTACCTTAAAGTCTGG
L	unitig1	+	unitig2	-	10M
P	contig1	unitig1+,unitig2-	10M
```

**References**: [GFA Format](https://gfa-spec.github.io/GFA-spec/GFA1.html)

## FASTG

**Extension**: `.fastg`

FASTG represents possible assembly paths using a FASTA-like structure that encodes graph topology.

Each node (`NODE_X`) in the graph (contig or unitig) contains:

* A path label (e.g. `:neighbor1,neighbor2...`) describing outgoing edges.
* A sequence entry representing the nucleotide sequence.

This is followed by a reverse-complement entry (`NODE_X'`). See the example below.

**Example**

```text
>NODE_1:NODE_2';
ATGCGTACGTTAG
>NODE_1';
CTAACGTACGCAT
>NODE_2:NODE_1',NODE_3';
CGGTAACCTGACC
>NODE_2';
GGTCAGGTTACCG
>NODE_3:NODE_2';
TTGACCGGATCGA
>NODE_3';
TCGATCCGGTCAA
```

**References**: [FASTG Format](https://web.archive.org/web/20211209213905/http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf)

**Note**: The FASTG format described in the official FASTG specification differs substantially from the FASTG files produced by MEGAHIT and early versions of SPAdes (prior to version 3.10.0).

## ASQG

**Extension**: `.asqg`

The ASQG format represents assembly graphs using overlapping contigs.

* Header records (start with the tag `HT`) denote metadata.
* Vertex records (start with the tag `VT`) denote the sequences.
* Edge description records (start with the tag `ED`) denote pairs of overlapping sequences.

**Example**

```text
HT	VN:i:1	ER:f:0	OL:i:45	IN:Z:reads.fa	CN:i:1	TE:i:0
VT	read1	GATCGATCTAGCTAGCTAGCTAGCTAGTTAGATGCATGCATGCTAGCTGG
VT	read2	CGATCTAGCTAGCTAGCTAGCTAGTTAGATGCATGCATGCTAGCTGGATA
VT	read3	ATCTAGCTAGCTAGCTAGCTAGTTAGATGCATGCATGCTAGCTGGATATT
ED	read2 read1 0 46 50 3 49 50 0 0
ED	read3 read2 0 47 50 2 49 50 0 0
```

**References**: [ASQG Format](https://github.com/jts/sga/wiki/ASQG-Format)

## DOT

**Extension**: `.dot` (GraphViz) or `.gv` (ABySS)

DOT files describe general graphs and can be used to render visual diagrams of assembly graphs.

**Example**

```text
graph {
  0 [
    id=0
    name=1
  ];
  1 [
    id=1
    name=2
  ];
  2 [
    id=2
    name=3
  ];
  3 [
    id=3
    name=4
  ];
  4 [
    id=4
    name=5
  ];

  1 -- 0;
  2 -- 1;
  3 -- 1;
  4 -- 2;
  4 -- 3;
}
```

To view it using GraphViz:

```bash
dot -Tpng graph.dot -o graph.png
```

**References**: [GraphViz DOT Format](https://graphviz.org/doc/info/lang.html)


The same example graph can be represented in ABySS DOT format, where:

* `l` denotes sequence length
* `d` denotes overlap length

**Example**

```text
digraph g {
"1+" [l=8]
"1-" [l=8]
"2+" [l=8]
"2-" [l=8]
"3+" [l=8]
"3-" [l=8]
"4+" [l=8]
"4-" [l=8]
"5+" [l=8]
"5-" [l=8]
"1+" -> "2+" [d=-7]
"2-" -> "1-" [d=-7]
"2+" -> "3+" [d=-7]
"3-" -> "2-" [d=-7]
"2+" -> "4+" [d=-7]
"4-" -> "2-" [d=-7]
"3+" -> "5+" [d=-7]
"5-" -> "3-" [d=-7]
"4+" -> "5+" [d=-7]
"5-" -> "4-" [d=-7]
}
```

**References**: [ABySS DOT Format](https://github.com/bcgsc/abyss/wiki/ABySS-File-Formats#dot)

---

## Comparison of Assembly Graph File Formats

| Feature | **GFA** | **FASTG** | **ASQG** | **DOT (GraphViz / ABySS)** |
|-------|--------|-----------|----------|-----------------------------|
| **Primary use** | General-purpose assembly graph representation | Compact representation of assembly paths | Read overlap graphs | Visualisation and graph layout |
| **Typical extension** | `.gfa` | `.fastg` | `.asqg` | `.dot`, `.gv` |
| **Graph type** | Directed graph | Directed graph | Directed graph | Directed or undirected |
| **Nodes represent** | Segments / unitigs | Contigs or unitigs | Reads or contigs | Generic graph nodes |
| **Edges represent** | Overlaps / links | Traversable paths | Sequence overlaps | Graph connections |
| **Stores sequences** | ✅ Yes | ✅ Yes | ✅ Yes | ❌ No |
| **Stores overlaps** | ✅ Yes | ✅ Implicitly | ✅ Yes | Optional |
| **Supports paths** | ✅ Yes (`P` records) | ✅ Yes (via node labels) | ❌ No | ❌ No |
| **Orientation-aware** | ✅ Yes | ✅ Yes | ✅ Yes | Optional |
| **Standardised** | ✅ Yes (GFA1) | ⚠️ Semi-standard | ❌ Tool-specific | ✅ GraphViz |
| **Used by assemblers** | SPAdes, Flye, hifiasm | MEGAHIT | SGA | ABySS |
| **Supported by agtools** | ✅ Full support | ✅ Supported | ✅ Supported | ✅ Supported |
| **Best for** | General graph analysis | Path-centric graphs | Overlap-based assemblies | Visualisation |
| **Visualisation support** | Via Bandage | Via Bandage | Limited | Native (GraphViz) |

---

## Summary

| Format | When to Use |
|------|-------------|
| **GFA** | Best all-purpose format for assembly graphs; recommended default |
| **FASTG** | Useful when working with assembler-generated paths |
| **ASQG** | Suitable for overlap-based assemblies |
| **DOT** | Best for visualisation and debugging |
