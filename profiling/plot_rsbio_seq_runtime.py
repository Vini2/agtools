"""Running time of loading a SPAdes assembly graph: biopython vs rsbio-seq.

"Loading" here is spades.get_contig_graph, which reads the GFA together with
the contigs FASTA and returns a ContigGraph.

Input: time_logs/rsbio_seq_spades_benchmark_res.csv, written by
rsbio_seq_benchmark.py (run once with `old` and once with `new` per dataset).
"""

import matplotlib.pyplot as plt
import pandas as pd

OLD_COLOUR, NEW_COLOUR = "orange", "blue"

df = pd.read_csv("time_logs/rsbio_seq_spades_benchmark_res.csv")
old = df[df["impl"] == "old"].sort_values("n_contigs")
new = df[df["impl"] == "new"].sort_values("n_contigs")

plt.figure(figsize=(5, 5))
plt.plot(old["n_contigs"] / 1e6, old["get_contig_graph_s"], "o-",
         color=OLD_COLOUR, label="biopython")
plt.plot(new["n_contigs"] / 1e6, new["get_contig_graph_s"], "o-",
         color=NEW_COLOUR, label="rsbio-seq")

plt.xlabel("Number of contigs (millions)")
plt.ylabel("Running time (s)")
plt.title("Running time - loading SPAdes assembly graph")
plt.grid(True)
plt.legend()
plt.savefig("plots/rsbio_seq_runtime_spades.pdf", dpi=600, bbox_inches="tight")
plt.savefig("plots/rsbio_seq_runtime_spades.png", dpi=150, bbox_inches="tight")
print("wrote plots/rsbio_seq_runtime_spades.{pdf,png}")
