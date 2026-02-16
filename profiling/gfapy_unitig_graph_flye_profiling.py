import cProfile
import os
import pstats
import time
import gc
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.ticker import FuncFormatter
from statistics import mean, stdev
from gfapy import Gfa
from pympler import asizeof  

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
    

def profile_call(func, *args, **kwargs):
    """
    Run func(*args, **kwargs) once under cProfile and return timings + deep memory.
    Memory is computed on the returned object via Pympler (recursive size).
    """
    gc.collect()  # reduce cross-run noise

    pr = cProfile.Profile()
    start_wall = time.perf_counter()
    pr.enable()
    obj = func(*args, **kwargs)  # capture object so we can size it
    pr.disable()
    end_wall = time.perf_counter()

    # Extract total cumulative time from profile stats
    ps = pstats.Stats(pr)
    cumtime = sum(stat[3] for stat in ps.stats.values())  # inclusive CPU time

    # Deep recursive size (bytes) of the returned object
    mem_bytes = asizeof.asizeof(obj)
    mem_mb = mem_bytes / (1024 * 1024)

    # Drop reference quickly to limit carryover between runs
    del obj
    gc.collect()

    return {
        "wall_time": end_wall - start_wall,
        "cumtime": cumtime,
        "mem_bytes": mem_bytes,
        "mem_MB": mem_mb,
        "profile": pr
    }

def profile_files(graph_files, runs=10):
    """Profile Gfa.from_file on each file multiple times."""
    results = []

    for file in graph_files:

        print(f"Profiling {file}...")
        
        wall_times = []
        cum_times = []
        mem_mbs = []  

        for _ in range(runs):
            res = profile_call(Gfa.from_file, file)
            wall_times.append(res["wall_time"])
            cum_times.append(res["cumtime"])
            mem_mbs.append(res["mem_MB"])  

        # GFA line stats
        count_S = grep_count("S", file)
        count_L = grep_count("L", file)
        count_P = grep_count("P", file)

        # File sizes in MB
        size_graph = get_file_size(file)

        results.append({
            "graph_file": file,
            # Wall clock
            "wall_min": min(wall_times),
            "wall_max": max(wall_times),
            "wall_mean": mean(wall_times),
            "wall_std": stdev(wall_times) if runs > 1 else 0.0,
            # CPU cumulative (from cProfile)
            "cum_min": min(cum_times),
            "cum_max": max(cum_times),
            "cum_mean": mean(cum_times),
            "cum_std": stdev(cum_times) if runs > 1 else 0.0,
            # Deep memory (MB) of returned object
            "mem_min_MB": min(mem_mbs),
            "mem_max_MB": max(mem_mbs),
            "mem_mean_MB": mean(mem_mbs),
            "mem_std_MB": stdev(mem_mbs) if runs > 1 else 0.0,
            # Context
            "gfa_S": count_S,
            "gfa_L": count_L,
            "gfa_P": count_P,
            "size_graph_MB": size_graph,
        })

    return pd.DataFrame(results)


def main():
    
    graph_files = [
        "/scratch/user/mall0133/Human_feces/SRR18490951/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/Human_feces/SRR18490961/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/Human_feces/SRR18491036/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/Human_feces/SRR18491148/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/Human_feces/SRR18491176/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/Human_feces/SRR18491204/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/Human_feces/SRR18491300/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/Human_feces/SRR18491309/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/Human_feces/SRR18491312/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/Human_feces/SRR18491319/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/MFD-LR/ERR10750395/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/MFD-LR/ERR11523645/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/MFD-LR/ERR11561019/flye_output/assembly_graph.gfa",
        "/scratch/user/mall0133/MFD-LR/ERR10750395/flye_output/assembly_graph.gfa"
    ]
    
    df = profile_files(graph_files, runs=10)

    # Save results
    df.to_csv("gfapy_unitig_graph_flye_profiling_res.csv", index=False)


if __name__ == "__main__":
    main()
