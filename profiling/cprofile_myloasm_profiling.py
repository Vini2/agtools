import cProfile
import os
import pstats
import time
from statistics import mean, stdev

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter

from agtools.assemblers import myloasm

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]


# GFA Line Count
def grep_count(line_prefix, file_path):
    if not os.path.exists(file_path):
        return 0
    with open(file_path, "r") as f:
        return sum(1 for line in f if line.startswith(line_prefix))


# Get file size in MB
def get_file_size(path):
    if os.path.exists(path):
        return os.path.getsize(path) / (1024 * 1024)  # Convert to MB
    return 0.0


def profile_call(func, *args, **kwargs):
    """Run func(*args, **kwargs) once under cProfile and return timings."""
    pr = cProfile.Profile()

    start_wall = time.perf_counter()
    pr.enable()
    func(*args, **kwargs)
    pr.disable()
    end_wall = time.perf_counter()

    # Extract total cumulative time from profile stats
    ps = pstats.Stats(pr)
    cumtime = sum(stat[3] for stat in ps.stats.values())  # inclusive CPU time

    return {"wall_time": end_wall - start_wall, "cumtime": cumtime, "profile": pr}


def profile_files(folder_paths, runs=10):
    """Profile myloasm.get_contig_graph on each folder multiple times."""
    results = []

    for folder_path in folder_paths:

        print(f"Profiling {folder_path}...")

        graph_file = os.path.join(folder_path, "final_contig_graph.gfa")
        contigs_file = os.path.join(folder_path, "assembly_primary.fa")

        wall_times = []
        cum_times = []

        for _ in range(runs):
            res = profile_call(myloasm.get_contig_graph, graph_file, contigs_file)
            wall_times.append(res["wall_time"])
            cum_times.append(res["cumtime"])

        # GFA line stats
        count_S = grep_count("S", graph_file)
        count_L = grep_count("L", graph_file)
        count_P = grep_count("P", graph_file)

        # File sizes in MB
        size_graph = get_file_size(graph_file)
        size_contigs = get_file_size(contigs_file)

        results.append(
            {
                "graph_file": folder_path,
                "wall_min": min(wall_times),
                "wall_max": max(wall_times),
                "wall_mean": mean(wall_times),
                "wall_std": stdev(wall_times) if runs > 1 else 0.0,
                "cum_min": min(cum_times),
                "cum_max": max(cum_times),
                "cum_mean": mean(cum_times),
                "cum_std": stdev(cum_times) if runs > 1 else 0.0,
                "gfa_S": count_S,
                "gfa_L": count_L,
                "gfa_P": count_P,
                "size_graph_MB": size_graph,
                "size_contigs_MB": size_contigs,
            }
        )

    return pd.DataFrame(results)


def main():

    folders = [
        "data/myloasm/ERR10750395",
        "data/myloasm/ERR11523645",
        "data/myloasm/ERR11561019",
        "data/myloasm/ERR11593880",
        "data/myloasm/ERR12040030",
    ]

    df = profile_files(folders, runs=10)

    # Save results
    df.to_csv("data/myloasm_cprofile_res.csv", index=False)

    # Read the profiling results from CSV file
    df_results = pd.read_csv("data/myloasm_cprofile_res.csv")

    # Plot running time with error bars
    # -----------------------------------------------------------
    x = (
        df_results["gfa_S"].to_numpy()
        + df_results["gfa_L"].to_numpy()
        + df_results["gfa_P"].to_numpy()
    )
    x_million = x / 1e6  # convert to millions
    y = df_results["wall_mean"].to_numpy()
    yerr = df_results["wall_std"].to_numpy()

    # Force regression through (0,0)
    m = (x_million * y).sum() / (x_million**2).sum()
    b = 0
    trend_y = m * x_million  # no intercept

    plt.figure(figsize=(8, 5))
    plt.errorbar(
        x_million,
        y,
        yerr=yerr,
        fmt="o",
        color="blue",
        ecolor="lightblue",
        elinewidth=2,
        capsize=4,
        label="Running time Mean ± Std",
    )

    # Plot trend line
    plt.plot(
        x_million, trend_y, "--", color="black", label=f"Trend line: y={m:.3f}x+{b:.3f}"
    )

    # Format x-axis with 1 decimal and "M"
    plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.1f}"))

    plt.xlabel("Number of GFA lines (millions)")
    plt.ylabel("Running time (s)")
    plt.title("Running time vs number of GFA lines for myloasm contig graph")
    plt.grid(True)

    # Save to file
    plt.savefig("plots/myloasm_time.png", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
