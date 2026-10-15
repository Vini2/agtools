import gc, os, resource, sys
impl, dataset = sys.argv[1], sys.argv[2]
sys.path.insert(0, f"{os.environ['HOME']}/{'baseline/src' if impl=='old' else 'mnt/agtools/src'}")
from agtools.assemblers import spades
from agtools.core.fasta_parser import FastaParser
D = f"{os.environ['HOME']}/mnt/agtools/profiling/data/SPAdes/{dataset}"
rss = lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
p = FastaParser(f"{D}/contigs.fasta"); gc.collect(); parser_only = rss()
del p; gc.collect()
g = spades.get_contig_graph(f"{D}/assembly_graph_after_simplification.gfa", f"{D}/contigs.fasta", f"{D}/contigs.paths")
gc.collect()
print(f"{impl},{dataset},{parser_only:.1f},{rss():.1f}")
