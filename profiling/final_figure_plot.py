import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter


configs = [
    {
        "name": "spades",
        "ag_file": "time_logs/agtools_unitig_graph_spades_profiling_res.csv",
        "gf_file": "time_logs/gfapy_unitig_graph_spades_profiling_res.csv",
        "runtime_title": "Running time - loading SPAdes GFA",
        "memory_title": "Memory usage - loading SPAdes GFA",
    },
    {
        "name": "megahit",
        "ag_file": "time_logs/agtools_unitig_graph_megahit_profiling_res.csv",
        "gf_file": "time_logs/gfapy_unitig_graph_megahit_profiling_res.csv",
        "runtime_title": "Running time - loading MEGAHIT converted GFA",
        "memory_title": "Memory usage - loading MEGAHIT converted GFA",
    },
    {
        "name": "flye",
        "ag_file": "time_logs/agtools_unitig_graph_flye_profiling_res.csv",
        "gf_file": "time_logs/gfapy_unitig_graph_flye_profiling_res.csv",
        "runtime_title": "Running time - loading Flye GFA",
        "memory_title": "Memory usage - loading Flye GFA",
    },
]


def plot_panel(ax, df_ag, df_gf, y_col, yerr_col, ylabel, title, panel_label):
    # AGTools
    x_ag = df_ag["gfa_S"] + df_ag["gfa_L"] + df_ag["gfa_P"]
    x_ag_m = x_ag / 1e6
    y_ag = df_ag[y_col]
    yerr_ag = df_ag[yerr_col]

    m_ag = (x_ag_m * y_ag).sum() / (x_ag_m**2).sum()
    ax.errorbar(
        x_ag_m,
        y_ag,
        yerr=yerr_ag,
        fmt="o",
        color="blue",
        ecolor="lightblue",
        capsize=4,
        label="agtools",
    )
    ax.plot(x_ag_m, m_ag * x_ag_m, "--", color="blue")

    # GfaPy
    x_gf = df_gf["gfa_S"] + df_gf["gfa_L"] + df_gf["gfa_P"]
    x_gf_m = x_gf / 1e6
    y_gf = df_gf[y_col]
    yerr_gf = df_gf[yerr_col]

    m_gf = (x_gf_m * y_gf).sum() / (x_gf_m**2).sum()
    ax.errorbar(
        x_gf_m,
        y_gf,
        yerr=yerr_gf,
        fmt="o",
        color="orange",
        ecolor="peachpuff",
        capsize=4,
        label="GfaPy",
    )
    ax.plot(x_gf_m, m_gf * x_gf_m, "--", color="orange")

    ax.set_xlabel("Number of GFA lines (millions)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.1f}"))
    ax.grid(True)
    ax.legend()

    # Panel label (A, B, C, ...)
    ax.text(
        -0.15,   # move left outside axes
        1.07,    # move above axes
        panel_label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
        ha="left",
    )


fig, axes = plt.subplots(2, 3, figsize=(15, 10))

panel_labels = ["A", "B", "C", "D", "E", "F"]
label_idx = 0

for col, cfg in enumerate(configs):
    df_ag = pd.read_csv(cfg["ag_file"])
    df_gf = pd.read_csv(cfg["gf_file"])

    plot_panel(
        axes[0, col],
        df_ag,
        df_gf,
        y_col="wall_mean",
        yerr_col="wall_std",
        ylabel="Running time (s)",
        title=cfg["runtime_title"],
        panel_label=panel_labels[label_idx],
    )
    label_idx += 1

    plot_panel(
        axes[1, col],
        df_ag,
        df_gf,
        y_col="mem_mean_MB",
        yerr_col="mem_std_MB",
        ylabel="Memory usage (MB)",
        title=cfg["memory_title"],
        panel_label=panel_labels[label_idx],
    )
    label_idx += 1


plt.tight_layout()
plt.savefig(
    "plots/runtime_memory_spades_megahit_flye_agtools_vs_gfapy.pdf",
    dpi=600,
    bbox_inches="tight",
)
plt.show()
