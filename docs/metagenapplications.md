# Using *agtools* in metagenomics

## Bin contigs by connected components

Contigs of a genome ideally form connected components in the assembly graph. Here is a minimal, straightforward example using *agtools* that “bins” contigs by connected component (each bin = one connected component in the contig graph). It uses the assembler helper to load a contig graph and then `get_connected_components()` to group contig.

```python
# 1) pick the assembler loader that matches your assembly
#    (SPAdes shown here; MEGAHIT/Flye/myloasm loaders work the same way)
from agtools.assemblers import spades
# If you used MEGAHIT:
#   from agtools.assemblers import megahit
# If you used Flye:
#   from agtools.assemblers import flye
# If you used myloasm:
#   from agtools.assemblers import myloasm

# --- files produced by your assembler ---
graph_file        = "assembly_graph_with_scaffolds.gfa"  # SPAdes GFA
contigs_fasta     = "contigs.fasta"                      # SPAdes contigs
contig_paths_file = "contigs.paths"                      # SPAdes contig paths

# 2) load the contig graph (ContigGraph)
cg = spades.get_contig_graph(graph_file, contigs_fasta, contig_paths_file)
# For MEGAHIT: cg = megahit.get_contig_graph(graph_file, contigs_fasta)
# For Flye:    cg = flye.get_contig_graph(graph_file, contigs_fasta, "assembly_info.txt")
# For myloasm: cg = myloasm.get_contig_graph(graph_file, contigs_fasta)

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
    This example bins purely by graph connectivity (topology). If you want coverage and nucleotide composition-aware binning, you have to combine these bins with additional heuristics or downstream methods. Please take a look at [MetaCoAG](https://doi.org/10.1007/978-3-031-04749-7_5).

## Identifying plasmid candidates

Here is a tiny, practical “plasmid candidate finder” using the *agtools* API. It flags circular contigs (self-loops) and reports their lengths—simple, fast, and a good first pass before deeper validation.

```python
# Choose the loader for your assembler (SPAdes shown here)
from agtools.assemblers import spades
from pathlib import Path

# --- files produced by your assembler ---
graph_file        = "assembly_graph_with_scaffolds.gfa"
contigs_fasta     = "contigs.fasta"
contig_paths_file = "contigs.paths"

# 1) Load contig graph
cg = spades.get_contig_graph(graph_file, contigs_fasta, contig_paths_file)

# 2) Grab circular contigs (self-loops) — these are strong plasmid candidates
circular_contigs = cg.self_loops  # list of contig names forming self-loops

# 3) Pull sequences, compute simple stats, and apply loose length filters
#    (plasmids are often a few kb to a few hundred kb; tweak as needed)
min_len = 1000
max_len = 500_000

candidates = []
for name in circular_contigs:
    seq = cg.get_contig_sequence(name)
    length = len(seq)
    if min_len <= length <= max_len:
        candidates.append({
            "name": name,
            "length_bp": length,
        })

# 4) Print a quick report
print(f"Found {len(candidates)} circular contig(s) in plausible plasmid size range [{min_len:,}-{max_len:,}] bp")
for i, rec in enumerate(sorted(candidates, key=lambda r: -r["length_bp"]), start=1):
    print(f"{i:2d}. {rec['name']}\tlen={rec['length_bp']:,} bp\tGC={rec['gc']:.3f}")

# 5) Optional: write results to a TSV
out = Path("plasmid_candidates.tsv")
with out.open("w") as fh:
    fh.write("contig_name\tlength_bp\tgc_fraction\n")
    for rec in candidates:
        fh.write(f"{rec['name']}\t{rec['length_bp']}\t{rec['gc']}\n")
print(f"\nSaved: {out.resolve()}")

```

!!! note
    Circular contigs do not means plasmids all the time. This example script is a first pass. For confirmation, you have to run gene/marker checks (replication proteins, MOB/relaxase, AMR markers) with downstream tools (e.g., PlasmidFinder or PLASMe) after you shortlist candidates from the graph.


## Identify bacteriophage candidates from an assembly graph

Bacteriophages (or phages) tend to form circular components in the assembly graph. Here is a simple “phage candidate finder” using *agtools*. It looks for simple cycles in the oriented unitig graph (i.e., circular paths) and keeps the ones whose estimated genome length falls in a typical bacteriophage range (default 10–300 kb; most sequenced phages cluster around 30–50 kb). You can adjust `MIN_LEN`/`MAX_LEN` as needed.

```python
# Strategy: find simple cycles in the oriented unitig graph and keep those in a phage-like size range.

# deps: pip install agtools python-igraph bidict
import igraph as ig
from bidict import bidict

from agtools.core.unitig_graph import UnitigGraph

# --- input: unitig GFA from your assembler (SPAdes/Flye/etc.) ---
gfa_path = "assembly_graph_with_scaffolds.gfa"

# size window for candidate phage genomes (kbp). Adjust as needed.
MIN_LEN = 10_000
MAX_LEN = 300_000

# 1) Load the unitig graph
ug = UnitigGraph.from_gfa(gfa_path)  # parses GFA; sequences fetched lazily from file

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
```