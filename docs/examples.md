# Usage examples

## Obtaining stats of an assembly graph

*agtools* can generate basic graph and sequence-based statistics for an assembly graph. You can use the `stats` subcommand provided through the command-line interface. Please refer to the [CLI reference](cli.md) for further details on the `stats` subcommand.

```bash
agtools stats -g tests/data/ESC/assembly_graph_with_scaffolds.gfa -o ./
```

A text file with the following details will be generated.

```text
Basic graph statistics for tests/data/ESC/assembly_graph_with_scaffolds.gfa:
Number of segments: 982
Number of links: 1265
Number of self-loops: 1
Number of connected components: 4
Average node degree: 2

Sequence-based statistics for tests/data/ESC/assembly_graph_with_scaffolds.gfa:
Total length of segments: 8337494 bp
Average segment length: 8490 bp
N50: 60706 bp
L50: 36 segment(s)
GC content: 43.52%
```

## Renaming elements of an assembly graph

*agtools* can rename segment IDs, path IDs and walk IDs by prepending a given prefix. You can use the `rename` subcommand provided through the command-line interface. Please refer to the [CLI reference](cli.md) for further details on the `rename` subcommand.

Here is an example GFA file.

```text
# GFA file: 10 seqments (20 bp) and 1 short seqment (10 bp: seqX)
H	VN:Z:1.0
# Segments
S	seq1	ATGCGTATGCGTATGCGTAA
S	seq2	CGTACTGACTGACTGACTGA
S	seq3	TGCATGCATGCATGCATGCA
S	seq4	ACGTACGTACGTACGTACGT
S	seq5	GTACGTACGTACGTACGTAC
S	seq6	CGTACGTACGTACGTACGTA
S	seq7	TACGTACGTACGTACGTACG
S	seq8	ATCGATCGATCGATCGATCG
S	seq9	GCGTGCGTGCGTGCGTGCGT
S	seq10	TTGCTTGCTTGCTTGCTTGC
S	seqX	ACGTACGTAC
L	seqX	+	seq2	-	5M
L	seq4	+	seq5	+	10M
L	seq6	-	seq7	+	7M
L	seq9	+	seq10	-	5M
J	seq5	+	seqX	-	*
J	seq3	-	seq8	+	*
J	seq7	+	seq2	+	*
C	seq1	+	seqX	-	10M	5	0
C	seq2	+	seqX	+	10M	2	0
P	seqpath1	seq1+,seqX+,seq3-	20M,10M,20M
P	seqpath2	seq4+,seq5+,seq6+	20M,20M,20M
W	seqread1	0	*	<seq6<seqX>seq9
W	seqread2	0	*	<seq1>seq2<seq3
```

You can run the following command to rename using the prefix `TEST`.

```bash
agtools rename -g tests/data/test_graph_1.gfa -p TEST -o ./
```

The renamed graph will look as follows.

```text
# GFA file: 10 seqments (20 bp) and 1 short seqment (10 bp: seqX)
H	VN:Z:1.0
# Segments
S	TEST_seq1	ATGCGTATGCGTATGCGTAA
S	TEST_seq2	CGTACTGACTGACTGACTGA
S	TEST_seq3	TGCATGCATGCATGCATGCA
S	TEST_seq4	ACGTACGTACGTACGTACGT
S	TEST_seq5	GTACGTACGTACGTACGTAC
S	TEST_seq6	CGTACGTACGTACGTACGTA
S	TEST_seq7	TACGTACGTACGTACGTACG
S	TEST_seq8	ATCGATCGATCGATCGATCG
S	TEST_seq9	GCGTGCGTGCGTGCGTGCGT
S	TEST_seq10	TTGCTTGCTTGCTTGCTTGC
S	TEST_seqX	ACGTACGTAC
L	TEST_seqX	+	TEST_seq2	-	5M
L	TEST_seq4	+	TEST_seq5	+	10M
L	TEST_seq6	-	TEST_seq7	+	7M
L	TEST_seq9	+	TEST_seq10	-	5M
J	TEST_seq5	+	TEST_seqX	-	*
J	TEST_seq3	-	TEST_seq8	+	*
J	TEST_seq7	+	TEST_seq2	+	*
C	TEST_seq1	+	TEST_seqX	-	10M	5	0
C	TEST_seq2	+	TEST_seqX	+	10M	2	0
P	TEST_seqpath1	TEST_seq1+,TEST_seqX+,TEST_seq3-	20M,10M,20M
P	TEST_seqpath2	TEST_seq4+,TEST_seq5+,TEST_seq6+	20M,20M,20M
W	TEST_seqread1	0	*	<TEST_seq6<TEST_seqX>TEST_seq9
W	TEST_seqread2	0	*	<TEST_seq1>TEST_seq2<TEST_seq3
```