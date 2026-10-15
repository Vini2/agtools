"""Peak memory when loading a SPAdes assembly graph: biopython vs rsbio-seq.

"Loading" here is spades.get_contig_graph, which reads the GFA together with
the contigs FASTA and returns a ContigGraph.

Input: time_logs/rsbio_seq_spades_memory_res.csv, written by
rsbio_seq_memory_probe.py (run once with `old` and once with `new` per
dataset).
"""

import matplotlib.pyplot as plt
import pandas as pd

OLD_COLOUR, NEW_COLOUR = "orange", "blue"

df = pd.read_csv("time_logs/rsbio_seq_spades_memory_res.csv")
old = df[df["impl"] == "old"].sort_values("n_contigs")
new = df[df["impl"] == "new"].sort_values("n_contigs")

plt.figure(figsize=(5, 5))
plt.plot(old["n_contigs"] / 1e6, old["after_graph_mb"], "o-",
         color=OLD_COLOUR, label="biopython")
plt.plot(new["n_contigs"] / 1e6, new["after_graph_mb"], "o-",
         color=NEW_COLOUR, label="rsbio-seq")

plt.xlabel("Number of contigs (millions)")
plt.ylabel("Peak memory (MB)")
plt.title("Peak memory - loading SPAdes assembly graph")
plt.ylim(bottom=0)
plt.grid(True)
plt.legend()
plt.savefig("plots/rsbio_seq_memory_spades.pdf", dpi=600, bbox_inches="tight")
plt.savefig("plots/rsbio_seq_memory_spades.png", dpi=150, bbox_inches="tight")
print("wrote plots/rsbio_seq_memory_spades.{pdf,png}")
