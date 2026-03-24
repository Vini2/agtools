## Concatenating multiple assembly graphs

*agtools* can concatenate multiple GFA files. This can be useful for downstream analysis of assemblies with multiple samples. You can use the `concat` subcommand provided through the command-line interface. Please refer to the [CLI reference](../cli.md) for further details on the `concat` subcommand.

!!! note
    The `concat` subcommand will NOT merge similar segments. It simply concatenates the segments and other elements in the GFA files.

You can run the following command to concatenate multiple GFA files. Make sure to provide each GFA file using `-g`.

```bash
agtools concat -g test_graph_1.gfa -g test_graph_2.gfa -g test_graph_3.gfa -o results/concatenated_graph.gfa
```

The lines with different tags `H`, `S`, `L`, `J`, `C`, `P` and `W` will be grouped together from the GFA files and written together.

!!! tip 
    If two GFA files have the same segment ID, path ID or walk ID, the run will fail. If you have duplicate IDs, please make sure to run the `rename` subcommand before concatenating.
