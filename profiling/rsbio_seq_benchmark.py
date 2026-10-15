import csv, os, pathlib, random, resource, sys, time

impl, dataset = sys.argv[1], sys.argv[2]
sys.path.insert(0, f"{os.environ['HOME']}/{'baseline/src' if impl == 'old' else 'mnt/agtools/src'}")

from agtools.assemblers import spades
from agtools.core.fasta_parser import FastaParser

DATA = pathlib.Path(os.environ["HOME"]) / "mnt/agtools/profiling/data/SPAdes" / dataset
graph = str(DATA / "assembly_graph_after_simplification.gfa")
contigs = str(DATA / "contigs.fasta")
paths = str(DATA / "contigs.paths")
size_mb = os.path.getsize(contigs) / 1048576

for suffix in (".fai", ".gzi"):
    try:
        os.remove(contigs + suffix)
    except FileNotFoundError:
        pass

def timed(fn):
    t = time.perf_counter(); value = fn(); return value, time.perf_counter() - t

parser, open_cold = timed(lambda: FastaParser(contigs))
_, open_warm = timed(lambda: FastaParser(contigs))

names = list(parser.index)
random.seed(0)
picks = [random.choice(names) for _ in range(1000)]
# First sequence access loads the index on the Rust side; time it separately
# from steady-state lookups.
_, reader_open = timed(lambda: getattr(parser, "reader", None))
bases, lookup = timed(lambda: sum(len(parser.get_sequence(n)) for n in picks))

rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
g, build = timed(lambda: spades.get_contig_graph(graph, contigs, paths))
rss_graph = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

if os.environ.get("STATS", "1") == "1":
    total_len, t_total = timed(g.calculate_total_length)
    n50, t_n50 = timed(g.calculate_n50_l50)
    gc, t_gc = timed(g.get_gc_content)
else:
    total_len, t_total = -1, float("nan")
    n50, t_n50 = (-1, -1), float("nan")
    gc, t_gc = -1.0, float("nan")

row = {
    "impl": impl, "dataset": dataset, "contigs_mb": round(size_mb, 1),
    "n_contigs": len(names), "vcount": g.vcount,
    "open_cold_s": round(open_cold, 3), "open_warm_s": round(open_warm, 3),
    "lookup_1000_s": round(lookup, 3), "reader_open_s": round(reader_open, 3), "get_contig_graph_s": round(build, 3),
    "total_length_s": round(t_total, 3), "n50_l50_s": round(t_n50, 3),
    "gc_content_s": round(t_gc, 3),
    "peak_rss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024),
    "rss_after_graph_mb": round(rss_graph),
    "bases_checksum": bases, "total_length": total_len,
    "n50": n50[0], "l50": n50[1], "gc": round(gc, 9),
}
out = pathlib.Path(os.environ["HOME"]) / "spades_bench3.csv"
write_header = not out.exists()
with out.open("a", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(row))
    if write_header:
        writer.writeheader()
    writer.writerow(row)
print(row)
