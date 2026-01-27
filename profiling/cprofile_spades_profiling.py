import cProfile
import os
import pstats
import time
from statistics import mean, stdev

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter

from agtools.assemblers import spades

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
    """Profile spades.get_contig_graph on each folder multiple times."""
    results = []

    for folder_path in folder_paths:

        print(f"Profiling {folder_path}...")

        graph_file = os.path.join(
            folder_path, "assembly_graph_after_simplification.gfa"
        )
        contigs_file = os.path.join(folder_path, "contigs.fasta")
        contig_paths_file = os.path.join(folder_path, "contigs.paths")

        wall_times = []
        cum_times = []

        for _ in range(runs):
            res = profile_call(
                spades.get_contig_graph, graph_file, contigs_file, contig_paths_file
            )
            wall_times.append(res["wall_time"])
            cum_times.append(res["cumtime"])

        # GFA line stats
        count_S = grep_count("S", graph_file)
        count_L = grep_count("L", graph_file)
        count_P = grep_count("P", graph_file)

        # File sizes in MB
        size_graph = get_file_size(graph_file)

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
            }
        )

    return pd.DataFrame(results)


def main():

    folders = [
        "data/SPAdes/ERR2752151",
        "data/SPAdes/ERR2752163",
        "data/SPAdes/ERR2752162",
        "data/SPAdes/ERR2752160",
        "data/SPAdes/ERR2752150",
        "data/SPAdes/ERR2752154",
        "data/SPAdes/ERR2752153",
        "data/SPAdes/ERR594355",
        "data/SPAdes/ERR599362",
        "data/SPAdes/ERR594375",
        "data/SPAdes/ERR594362",
        "data/SPAdes/ERR594361",
        "data/SPAdes/ERR2750828",
        "data/SPAdes/ERR599357",
        "data/SPAdes/ERR2750826",
        "data/SPAdes/ERR594360",
        "data/SPAdes/ERR599383",
        "data/SPAdes/ERR599370",
        "data/SPAdes/ERR321017",
        "data/SPAdes/ERR321018",
    ]

    df = profile_files(folders, runs=10)

    # Save results
    df.to_csv("data/spades_cprofile_res.csv", index=False)

    # Read the profiling results from CSV file
    df_results = pd.read_csv("data/spades_cprofile_res.csv")

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

    plt.figure(figsize=(6, 3))
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
    plt.title("SPAdes contig graph")
    plt.grid(True)

    # Save to file
    plt.savefig("plots/spades_time.png", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
