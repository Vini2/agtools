# Example applications of *agtools*

## Bin contigs by connected components

Contigs from a genome ideally form connected components in the assembly graph. Below is a minimal, straightforward example using *agtools* that “bins” contigs by connected components (each bin = one connected component in the contig graph). It uses *agtools* to load a contig graph and then obtain connected components through the `get_connected_components()` function to group contigs.

```python
# 1) pick the assembler loader that matches your assembly
#    (metaSPAdes shown here)
from agtools.assemblers import spades

# --- files produced by your assembler ---
graph_file        = "ERR321017_assembly_graph.gfa"  # metaSPAdes GFA
contigs_fasta     = "ERR321017_contigs.fasta"       # metaSPAdes contigs
contig_paths_file = "ERR321017_contigs.paths"       # metaSPAdes contig paths

# 2) load the contig graph (ContigGraph)
cg = spades.get_contig_graph(graph_file, contigs_fasta, contig_paths_file)

# 3) compute connected components (bins) at contig level
components = cg.get_connected_components()   # list of lists of internal IDs

# 4) map internal IDs -> contig names and collect bins
#    (contig_names[i] gives the name for the internal ID i)
bins = []
for comp in components:
    if len(comp) > 1: # ignore isolated contigs
        bin_names = [cg.contig_names[i] for i in comp]
        bins.append(bin_names)

# 5) report / save
print(f"Found {len(bins)} bins")
for k, bin_names in enumerate(bins, start=1):
    print(f"Bin {k}  (n={len(bin_names)}): first few -> {bin_names[:5]}")

# Optional: write one TSV (bin_id, contig_name)
with open("contig_bins.tsv", "w") as out:
    out.write("bin_id\tcontig_name\n")
    for k, bin_names in enumerate(bins, start=1):
        for name in bin_names:
            out.write(f"{k}\t{name}\n")

```

!!! note
    This example bins purely by graph connectivity (topology). If you want coverage- and nucleotide composition-aware binning, you will need to combine these bins with additional heuristics or downstream methods. Please refer to [MetaCoAG](https://doi.org/10.1007/978-3-031-04749-7_5).

<div class="button-container">
  <a href="../resources/ERR321017_assembly.zip" download class="download-button">
    <i class="bi bi-download"></i> Download the example input data: ERR321017_assembly.zip
  </a>
</div>

<hr>


## Identify plasmid candidates

Below is a simple example “plasmid candidate finder” using the *agtools* API. It flags circular contigs (self-loops) and reports their lengths. This is a useful first pass before deeper validation.

```python
from agtools.core.unitig_graph import UnitigGraph
from pathlib import Path

# --- input: unitig GFA from your assembler (SPAdes/Flye/etc.) ---
gfa_path = "ERR10750395_assembly_graph.gfa"

# 1) Load the unitig graph
ug = UnitigGraph.from_gfa(gfa_path)

# 2) Grab circular unitigs (self-loops) — these are strong plasmid candidates
circular_unitigs = ug.self_loops  # list of unitig names forming self-loops

# 3) Pull sequences, compute simple stats, and apply loose length filters
#    (plasmids are often a few kb to a few hundred kb; tweak as needed)
min_len = 1000
max_len = 500_000

candidates = []
for name in circular_unitigs:
    seq = ug.get_segment_sequence(ug.segment_names[name])
    length = len(seq)
    if min_len <= length <= max_len:
        candidates.append({
            "name": name,
            "length_bp": length,
        })

# 4) Print a quick report
print(f"Found {len(candidates)} circular unitig(s) in plausible plasmid size range [{min_len:,}-{max_len:,}] bp")
for i, rec in enumerate(sorted(candidates, key=lambda r: -r["length_bp"]), start=1):
    print(f"{i:2d}. {rec['name']}\tlen={rec['length_bp']:,} bp")

# 5) Optional: write results to a TSV
out = Path("plasmid_candidates.tsv")
with out.open("w") as fh:
    fh.write("unitig_name\tlength_bp\n")
    for rec in candidates:
        fh.write(f"{rec['name']}\t{rec['length_bp']}\n")
print(f"\nSaved: {out.resolve()}")

```

!!! note
    Circular segments do not always indicate plasmids. This script is a first pass. For confirmation, run gene or marker checks (e.g. replication proteins, MOB/relaxase, AMR markers) with downstream tools (such as PlasmidFinder or PLASMe) after shortlisting candidates from the graph.

<div class="button-container">
  <a href="../resources/ERR10750395_assembly.zip" download class="download-button">
    <i class="bi bi-download"></i> Download the example input data: ERR10750395_assembly.zip
  </a>
</div>

<hr>


## Identify bacteriophage candidates

Bacteriophages (or phages) with circular genomes often form circular components in the assembly graph. Below is a simple “phage candidate finder” using *agtools*. It looks for simple cycles in the oriented unitig graph (i.e. circular paths) and keeps the ones whose estimated genome length falls within a typical bacteriophage range (default 10–300 kb; many sequenced phages cluster around 30–50 kb). You can adjust `MIN_LEN` and `MAX_LEN` as needed.

```python
import igraph as ig
from bidict import bidict

from agtools.core.unitig_graph import UnitigGraph

# --- input: unitig GFA from your assembler (SPAdes/Flye/etc.) ---
gfa_path = "ERR10359653_component_graph.gfa"

# size window for candidate phage genomes (kbp). Adjust as needed.
MIN_LEN = 10_000
MAX_LEN = 300_000

# 1) Load the unitig graph
ug = UnitigGraph.from_gfa(gfa_path)

# 2) Build an oriented (directed) version of the graph, as in the agtools example
#    Each segment gets a forward (+) and reverse (-) node; edges carry orientation.
oriented_nodes = bidict()
i = 0
for seg_name in ug.segment_names:   # internalID -> name mapping
    oriented_nodes[f"{seg_name}+"] = i
    oriented_nodes[f"{seg_name}-"] = i + 1
    i += 2

directed_ug = ig.Graph(directed=True)
directed_ug.add_vertices(len(oriented_nodes))

# label vertices with their oriented names for easy reading
oriented_nodes_rev = oriented_nodes.inverse
for vid in range(len(oriented_nodes)):
    directed_ug.vs[vid]["name"] = oriented_nodes_rev[vid]

# Oriented edges come from ug.oriented_links
edge_list = []
for from_id in ug.oriented_links:
    for to_id in ug.oriented_links[from_id]:
        for (from_or, to_or) in ug.oriented_links[from_id][to_id]:
            source = f"{ug.segment_names[from_id]}{from_or}"
            target = f"{ug.segment_names[to_id]}{to_or}"
            edge_list.append((source, target))
directed_ug.add_edges([(oriented_nodes[s], oriented_nodes[t]) for s, t in edge_list])

# 3) Enumerate simple cycles (each is a circular path with no repeated vertices)
cycles = directed_ug.simple_cycles()

def oriented_name_to_parts(oname):
    # "9001+" -> ("9001", "+")
    return oname[:-1], oname[-1]

def cycle_length_bp(cycle_vertices):
    # Sum unique segment lengths, minus overlaps along each oriented edge in the cycle
    oriented_names = [directed_ug.vs[v]["name"] for v in cycle_vertices]
    segs = []
    for on in oriented_names:
        sname, _ = oriented_name_to_parts(on)
        segs.append(sname)
    unique_seg_sum = sum(ug.segment_lengths[s] for s in set(segs))

    # sum overlaps along edges (wrap around end->start)
    overlap_sum = 0
    for a, b in zip(oriented_names, oriented_names[1:] + oriented_names[:1]):
        sa, oa = oriented_name_to_parts(a)
        sb, ob = oriented_name_to_parts(b)
        key = (ug.segment_name_to_id[sa], oa, ug.segment_name_to_id[sb], ob)
        overlap_sum += ug.link_overlap.get(key, 0)

    return max(0, unique_seg_sum - overlap_sum)

# 4) Score, filter, and lightly deduplicate cycles
seen_signatures = set()
candidates = []
for cyc in cycles:
    # signature: set of (segment names) ignoring orientation, to collapse reverse-complement/rotations
    oriented_names = [directed_ug.vs[v]["name"] for v in cyc]
    base_names = tuple(sorted({on[:-1] for on in oriented_names}))
    if base_names in seen_signatures:
        continue
    seen_signatures.add(base_names)

    L = cycle_length_bp(cyc)
    if MIN_LEN <= L <= MAX_LEN:
        candidates.append({
            "length_bp": L,
            "n_segments": len(base_names),
            "path_oriented": " -> ".join(oriented_names)
        })

# 5) Report
candidates.sort(key=lambda r: -r["length_bp"])
print(f"Phage-like circular candidates: {len(candidates)} (size window {MIN_LEN:,}-{MAX_LEN:,} bp)")
for i, rec in enumerate(candidates, 1):
    print(f"{i:2d}. ~{rec['length_bp']:,} bp  | segments: {rec['n_segments']}\n    {rec['path_oriented']}\n")

# Optional: write results to a TSV
out = Path("phage_like_candidates.tsv")
with out.open("w") as fh:
    fh.write("path_number\tlength_bp\tn_segments\tpath_oriented\n")
    for i, rec in enumerate(candidates, 1):
        fh.write(f"{i}\t{rec['length_bp']}\t{rec['n_segments']}\t{rec['path_oriented']}\n")
print(f"\nSaved: {out.resolve()}")
```

!!! note
    This is a graph-only heuristic. Circular sequences do not always indicate phages (they could be plasmids, for example). For biological confirmation, follow up with hallmark phage gene searches (terminase large subunit, capsid, portal, tail) using tools such as [Pharokka](https://github.com/gbouras13/pharokka), [Phold](https://github.com/gbouras13/phold) and [Phynteny](https://github.com/susiegriggo/Phynteny_transformer) after shortlisting candidates. If you want to recover both circular and linear phage genomes more accurately, please refer to [Phables](https://github.com/Vini2/phables).

<div class="button-container">
  <a href="../resources/ERR10359653_assembly.zip" download class="download-button">
    <i class="bi bi-download"></i> Download the example input data: ERR10359653_assembly.zip
  </a>
</div>

<hr>


## Haplotype phasing from assembly graph bubbles

Haplotype phasing aims to distinguish and reconstruct the individual haplotypes present in a diploid or polyploid genome. In assembly graphs, this can often be approached by identifying and resolving **bubbles**, which are alternative paths in the graph corresponding to sequence differences between haplotypes.

Below is a minimal local haplotype phasing example using graph bubbles with the agtools API. It:

1. loads a `UnitigGraph` from a GFA,
2. builds the **oriented unitig graph** (the +/− version),
3. detects simple **bubbles** (two vertex-disjoint paths that split at a source and re-join at a sink), and
4. outputs the two **allele sequences** per bubble by stitching unitig sequences while subtracting per-edge overlaps.


```python
import igraph as ig
from bidict import bidict

from agtools.core.unitig_graph import UnitigGraph

# ---- Input: unitig GFA from your assembler (SPAdes/Flye/etc.) ----
gfa_path = "chr22_toy_bubbles.gfa"

# 1) Load unitig graph (agtools)
ug = UnitigGraph.from_gfa(gfa_path)

# 2) Build oriented (+/−) unitig graph
oriented_nodes = bidict()
vid = 0
for seg_name in ug.segment_names:
    oriented_nodes[f"{seg_name}+"] = vid
    oriented_nodes[f"{seg_name}-"] = vid + 1
    vid += 2

G = ig.Graph(directed=True)
G.add_vertices(len(oriented_nodes))

# label vertices with oriented names for convenience
revmap = oriented_nodes.inverse
for v in range(G.vcount()):
    G.vs[v]["name"] = revmap[v]

# add oriented edges using ug.oriented_links
edges = []
for from_id in ug.oriented_links:
    for to_id in ug.oriented_links[from_id]:
        for (from_or, to_or) in ug.oriented_links[from_id][to_id]:
            s = f"{ug.segment_names[from_id]}{from_or}"
            t = f"{ug.segment_names[to_id]}{to_or}"
            edges.append((oriented_nodes[s], oriented_nodes[t]))
G.add_edges(edges)

# 3) Utilities for stitching sequences along an oriented path

def oriented_name_parts(oname):
    return oname[:-1], oname[-1]  # ("9001", "+") etc.

def oriented_seq(oname):
    seg, orient = oriented_name_parts(oname)
    seq = ug.get_segment_sequence(seg)  # returns Bio.Seq.Seq
    return str(seq if orient == "+" else seq.reverse_complement())

def stitch_oriented_path(v_ids):
    """Concatenate sequences for an oriented vertex path, subtracting edge overlaps."""
    names = [G.vs[v]["name"] for v in v_ids]
    if not names:
        return ""
    out = oriented_seq(names[0])
    for a, b in zip(names, names[1:]):
        sa, oa = oriented_name_parts(a)
        sb, ob = oriented_name_parts(b)
        key = (ug.segment_name_to_id[sa], oa, ug.segment_name_to_id[sb], ob)
        ov = ug.link_overlap.get(key, 0)  # overlap in bp (agtools stores this)
        out += oriented_seq(b)[ov:]
    return out

# 4) Find simple "bubbles" and phase them locally
# Heuristic: for any vertex with out-degree >=2, pick two distinct successors,
# find the first reconvergence vertex t reachable from both, then take
# the shortest disjoint paths s->t (one through each successor).
MAX_DEPTH = 80     # limit search for reconvergence
MAX_BUBBLES = 10   # cap for a quick demo

def first_common_reconvergence(n1, n2, depth=MAX_DEPTH):
    # BFS distance dicts from each start
    d1 = G.shortest_paths_dijkstra(source=n1, target=None)[0]
    d2 = G.shortest_paths_dijkstra(source=n2, target=None)[0]
    candidates = [v for v in range(G.vcount())
                  if d1[v] < float("inf") and d2[v] < float("inf")
                  and d1[v] <= depth and d2[v] <= depth
                  and G.indegree(v) >= 2]
    if not candidates:
        return None
    # prefer earliest joint target by summed distance
    return min(candidates, key=lambda v: (d1[v] + d2[v]))

bubbles = []
for s in range(G.vcount()):
    succ = G.successors(s)
    if len(succ) < 2:
        continue
    # try pairs of outgoing branches
    for i in range(len(succ)):
        for j in range(i + 1, len(succ)):
            n1, n2 = succ[i], succ[j]
            t = first_common_reconvergence(n1, n2)
            if t is None:
                continue
            # shortest paths s->t that begin with the chosen successor
            p1 = [s] + G.get_shortest_paths(n1, to=t, output="vpath")[0]
            p2 = [s] + G.get_shortest_paths(n2, to=t, output="vpath")[0]
            # ensure paths meet only at s and t (vertex-disjoint interiors)
            if set(p1[1:-1]).isdisjoint(set(p2[1:-1])):
                bubbles.append((s, t, p1, p2))
            if len(bubbles) >= MAX_BUBBLES:
                break
        if len(bubbles) >= MAX_BUBBLES:
            break
    if len(bubbles) >= MAX_BUBBLES:
        break

# 5) Emit phased alleles (two sequences per bubble)
print(f"Found {len(bubbles)} bubble(s) suitable for local phasing")
for idx, (s, t, pA, pB) in enumerate(bubbles, 1):
    hapA = stitch_oriented_path(pA)
    hapB = stitch_oriented_path(pB)
    s_name = G.vs[s]["name"]; t_name = G.vs[t]["name"]
    print(f"\nBubble {idx}: {s_name}  →  {t_name}")
    print(f" allele A length: {len(hapA):,} bp")
    print(f" allele B length: {len(hapB):,} bp")

    # Optional: write FASTA for downstream polishing/validation
    # with open(f"bubble_{idx}_alleleA.fasta", "w") as fh:
    #     fh.write(f">bubble_{idx}_alleleA\n{hapA}\n")
    # with open(f"bubble_{idx}_alleleB.fasta", "w") as fh:
    #     fh.write(f">bubble_{idx}_alleleB\n{hapB}\n")

```

!!! note
    This is a minimal example intended to demonstrate the utility of *agtools*. For chromosome-scale phasing, you will need to incorporate read or linkage evidence (e.g. long reads, Hi-C) to chain many bubbles consistently.

<div class="button-container">
  <a href="../resources/chr22_toy_bubbles.gfa" download class="download-button">
    <i class="bi bi-download"></i> Download the example input data: chr22_toy_bubbles.gfa
  </a>
</div>

<hr>