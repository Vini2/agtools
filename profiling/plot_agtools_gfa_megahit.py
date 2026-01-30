import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# -----------------------------------------------------------
# Load data
# -----------------------------------------------------------
df_ag = pd.read_csv("time_logs/agtools_unitig_graph_megahit_profiling_res.csv")
df_gf = pd.read_csv("time_logs/gfapy_unitig_graph_megahit_profiling_res.csv")

# ===========================================================
# FIGURE 1 — RUNTIME
# ===========================================================
plt.figure(figsize=(6, 3))

# AGTools
x_ag = df_ag["gfa_S"] + df_ag["gfa_L"] + df_ag["gfa_P"]
x_ag_m = x_ag / 1e6
y_ag = df_ag["wall_mean"]
yerr_ag = df_ag["wall_std"]

m_ag = (x_ag_m * y_ag).sum() / (x_ag_m**2).sum()
plt.errorbar(
    x_ag_m, y_ag, yerr=yerr_ag,
    fmt="o", color="blue", ecolor="lightblue",
    capsize=4, label="agtools"
)
plt.plot(x_ag_m, m_ag * x_ag_m, "--", color="blue")

# Gfapy
x_gf = df_gf["gfa_S"] + df_gf["gfa_L"] + df_gf["gfa_P"]
x_gf_m = x_gf / 1e6
y_gf = df_gf["wall_mean"]
yerr_gf = df_gf["wall_std"]

m_gf = (x_gf_m * y_gf).sum() / (x_gf_m**2).sum()
plt.errorbar(
    x_gf_m, y_gf, yerr=yerr_gf,
    fmt="o", color="orange", ecolor="peachpuff",
    capsize=4, label="Gfapy"
)
plt.plot(x_gf_m, m_gf * x_gf_m, "--", color="orange")

plt.xlabel("Number of GFA lines (millions)")
plt.ylabel("Running time (s)")
plt.title("agtools vs. Gfapy MEGAHIT")
plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.1f}"))
plt.grid(True)
plt.legend()
plt.savefig("plots/runtime_megahit_agtools_vs_gfapy.png", dpi=300, bbox_inches="tight")
plt.show()


# ===========================================================
# FIGURE 2 — MEMORY
# ===========================================================
plt.figure(figsize=(6, 3))

# AGTools
y_ag_mem = df_ag["mem_mean_MB"]
yerr_ag_mem = df_ag["mem_std_MB"]

m_ag_mem = (x_ag_m * y_ag_mem).sum() / (x_ag_m**2).sum()
plt.errorbar(
    x_ag_m, y_ag_mem, yerr=yerr_ag_mem,
    fmt="o", color="blue", ecolor="lightblue",
    capsize=4, label="agtools"
)
plt.plot(x_ag_m, m_ag_mem * x_ag_m, "--", color="blue")

# Gfapy
y_gf_mem = df_gf["mem_mean_MB"]
yerr_gf_mem = df_gf["mem_std_MB"]

m_gf_mem = (x_gf_m * y_gf_mem).sum() / (x_gf_m**2).sum()
plt.errorbar(
    x_gf_m, y_gf_mem, yerr=yerr_gf_mem,
    fmt="o", color="orange", ecolor="peachpuff",
    capsize=4, label="Gfapy"
)
plt.plot(x_gf_m, m_gf_mem * x_gf_m, "--", color="orange")

plt.xlabel("Number of GFA lines (millions)")
plt.ylabel("Memory usage (MB)")
plt.title("agtools vs. Gfapy MEGAHIT")
plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.1f}"))
plt.grid(True)
plt.legend()
plt.savefig("plots/memory_megahit_agtools_vs_gfapy.png", dpi=300, bbox_inches="tight")
plt.show()
