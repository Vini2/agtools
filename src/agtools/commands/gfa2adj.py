#!/usr/bin/env python3

import pandas as pd

from agtools.core.graph import UnitigGraph
from agtools.log_config import logger


def gfa2adj(gfa_file, output_path):

    ug = UnitigGraph.from_gfa(gfa_file)

    adj_matrix = ug.graph.get_adjacency()

    labels = list(ug.segment_names.values())

    adj_df = pd.DataFrame(adj_matrix, index=labels, columns=labels)

    output_file = f"{output_path}/adjacency_matrix.tsv"
    adj_df.to_csv(output_file, sep="\t")

    return output_file


# TODO: show isolated segments in the adjacency matrix
