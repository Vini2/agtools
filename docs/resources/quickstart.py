from agtools.core.unitig_graph import UnitigGraph

import igraph as ig

# Load the unitig-level assembly graph from the GFA file
ug = UnitigGraph.from_gfa("assembly_graph.gfa")

# Print basic graph-based statistics
print(f"Number of vertices (segments): {ug.vcount}")
print(f"Number of edges (links): {ug.ecount}")
print(f"Number of self-loops (repeats): {ug.self_loops}")
print(f"Average node degree: {ug.calculate_average_node_degree()}")
print(f"Number of connected components: {len(ug.get_connected_components())}")

# Print basic sequence-based statistics
print(f"Total length of segments: {ug.calculate_total_length()}")
print(f"Average length of segments: {ug.calculate_average_segment_length()}")
print(f"GC content of segments: {ug.get_gc_content()}")
print(f"N50 and L50: {ug.calculate_n50_l50()}")

# Print oriented links
for from_link in ug.oriented_links:
    for to_link in ug.oriented_links[from_link]:
        for orient in ug.oriented_links[from_link][to_link]:
            print(from_link, orient[0], "->", to_link, orient[1])

# Print paths
for path in ug.paths:
    # Path name: path string, path overlaps
    print(f"{path}: {ug.paths[path][0]}\t{ug.paths[path][1]}")

# Get neighbours of a segment
print(f"Neighbours of segment 5: {ug.get_neighbors("5")}")

# Print adjacency matrix
print("Adjancency matrix:")
print(ug.get_adjacency_matrix())

# Print segments
for seg_id in ug.segment_names:
    # internal segment ID, segment name, segment sequence, segment length
    print(seg_id, ug.segment_names[seg_id], ug.get_segment_sequence(ug.segment_names[seg_id]), ug.segment_lengths[ug.segment_names[seg_id]])

# Plot the graph using igraph
ig.plot(
    ug.graph,                           # graph object
    "graph_plot.png",                   # file name
    vertex_label=ug.graph.vs["name"],   # label names
    vertex_size=40,                     # vertex size
    vertex_frame_width=2.0,             #vertex frame width
    vertex_label_size=20.0,             # vertex label size
)

