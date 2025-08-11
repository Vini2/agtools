#!/usr/bin/env python3

import gc
import os
import time
import statistics
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from memory_profiler import memory_usage
from agtools.assemblers import spades

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2025, agtools Project"
__credits__ = ["Vijini Mallawaarachchi"]

# GFA Line Count
def grep_count(line_prefix, file_path):
    if not os.path.exists(file_path):
        return 0
    with open(file_path, 'r') as f:
        return sum(1 for line in f if line.startswith(line_prefix))

# Get file size in MB
def get_file_size(path):
    if os.path.exists(path):
        return os.path.getsize(path) / (1024 * 1024)  # Convert to MB
    return 0.0

# Main profiling code
def profile_get_contig_graph(graph_file, contigs_file, contig_paths_file):
    start = time.perf_counter()
    
    mem_trace, _ = memory_usage(
        (spades.get_contig_graph, (graph_file, contigs_file, contig_paths_file)),
        retval=True,
        interval=0.01,
        timeout=None
    )
    
    end = time.perf_counter()
    
    elapsed = end - start
    peak_mem = max(mem_trace) - min(mem_trace)
    
    return elapsed, peak_mem

def profile_folder(folder_path, runs=10):
    graph_file = os.path.join(folder_path, "assembly_graph_after_simplification.gfa")
    contigs_file = os.path.join(folder_path, "contigs.fasta")
    contig_paths_file = os.path.join(folder_path, "contigs.paths")

    times = []
    peak_mems = []

    for _ in range(runs):
        elapsed, peak = profile_get_contig_graph(graph_file, contigs_file, contig_paths_file)
        times.append(elapsed)
        peak_mems.append(peak)

    # GFA line stats
    count_S = grep_count("S", graph_file)
    count_L = grep_count("L", graph_file)
    count_P = grep_count("P", graph_file)

    # File sizes in MB
    size_graph = get_file_size(graph_file)
    size_fasta = get_file_size(contigs_file)
    size_paths = get_file_size(contig_paths_file) 

    return {
        "folder": os.path.basename(folder_path),
        "time_min": min(times),
        "time_max": max(times),
        "time_mean": statistics.mean(times),
        "time_std": statistics.stdev(times),
        "peak_mem_min": min(peak_mems),
        "peak_mem_max": max(peak_mems),
        "peak_mem_mean": statistics.mean(peak_mems),
        "peak_mem_std": statistics.stdev(peak_mems),
        "gfa_S": count_S,
        "gfa_L": count_L,
        "gfa_P": count_P,
        "size_graph_MB": size_graph,
        "size_fasta_MB": size_fasta,
        "size_paths_MB": size_paths,
    }

def batch_profile(folders, runs=10):
    results = []
    
    for folder in folders:
        
        print(f"Profiling: {folder}")
        
        stats = profile_folder(folder, runs)
        results.append(stats)
        
        # Manually trigger garbage collection
        gc.collect()
        
    return pd.DataFrame(results)

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

# Run profiling
df_results = batch_profile(folders, runs=20)

# Save to CSV
df_results.to_csv("profiling_spades.csv", index=False)

# Read the profiling results from CSV file
df_results = pd.read_csv("data/profiling_spades.csv")

# Plot running time with error bars
# -----------------------------------------------------------
x = df_results["size_graph_MB"]
y = df_results["time_mean"]
yerr = df_results["time_std"]

# Fit a simple linear regression for the trend line
m, b = np.polyfit(x, y, 1)  # slope, intercept
trend_y = m * x + b

plt.figure(figsize=(8, 5))
plt.errorbar(x, y, yerr=yerr, fmt='o', color='blue',
             ecolor='lightblue', elinewidth=2, capsize=4,
             label='Running time Mean ± Std')

# Plot trend line
plt.plot(x, trend_y, '--', color='black', label=f'Trend line: y={m:.3f}x+{b:.3f}')

plt.xlabel("Size of the graph file (MB)")
plt.ylabel("Running time (s)")
plt.title("Running time vs size of the graph file for SPAdes contig graph")
plt.legend()
plt.grid(True)

# Save to file
plt.savefig("plots/spades_time.png", dpi=300, bbox_inches='tight')
plt.show()

# Plot Peak Memory with error bars
# -----------------------------------------------------------
x = df_results["size_graph_MB"]
y = df_results["peak_mem_mean"]
yerr = df_results["peak_mem_std"]

# Fit a simple linear regression for the trend line
m, b = np.polyfit(x, y, 1)  # slope, intercept
trend_y = m * x + b

plt.figure(figsize=(8, 5))
plt.errorbar(x, y, yerr=yerr, fmt='o', color='red',
             ecolor='lightblue', elinewidth=2, capsize=4,
             label='Peak Memory Mean ± Std')

# Plot trend line
plt.plot(x, trend_y, '--', color='black', label=f'Trend line: y={m:.3f}x+{b:.3f}')

plt.xlabel("Size of the graph file (MB)")
plt.ylabel("Peak Memory (MB)")
plt.title("Peak memory vs size of the graph file for SPAdes contig graph")
plt.legend()
plt.grid(True)

# Save to file
plt.savefig("plots/spades_mem.png", dpi=300, bbox_inches='tight')
plt.show()
