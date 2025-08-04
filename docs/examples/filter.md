## Filtering segments from an assembly graph

*agtools* can filter sequences given a minimum segment length. Sequences that are shorted than minimum length will be removed along with any other elements that contain these segments. You can use the `filter` subcommand provided through the command-line interface. Please refer to the [CLI reference](cli.md) for further details on the `filter` subcommand.

Here is an example GFA file.

```text
# Test GFA file
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

We want to remove segments shorted than 15bp. You can run the following command to remove segments shorted than 15bp.

```bash
agtools filter -g test_graph.gfa -l 15 -o ./
```

The filtered graph file will look as follows.

```text
# Test GFA file
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
L	seq4	+	seq5	+	10M
L	seq6	-	seq7	+	7M
L	seq9	+	seq10	-	5M
J	seq3	-	seq8	+	*
J	seq7	+	seq2	+	*
P	seqpath1	seq1+,seqX+,seq3-	20M,10M,20M
P	seqpath2	seq4+,seq5+,seq6+	20M,20M,20M
W	seqread2	0	*	<seq1>seq2<seq3
```

!!! note
    Note that `seqX` and all links, jumps, containments, paths and walks containing `seqX` have been removed from the assembly graph.