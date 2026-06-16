# Standardized precomputed community files

Each file in `data/community_detection/precomputed/` contains one variable:

```matlab
community_results
```

Important fields:

- `network_name`: standardized network name.
- `source_label`: label from the original analysis file.
- `gamma`: second-pass gamma values.
- `consensus_partitions`: nodes × gamma consensus partitions.
- `n_communities`: number of raw consensus communities at each gamma.
- `zrand_mean`, `zrand_variance`: stability metrics at each gamma.
- `zrand_valid`: valid real-valued stability entries.
- `modularity`: gamma × repetition modularity values.
- `mean_modularity`: mean modularity at each gamma.
- `input_matrix`: original edge-weight matrix used in the analysis.
- `transformed_matrix`: positive z-shifted matrix supplied to Louvain.
- `paper.selected_index`: one-based selected index in `gamma`.
- `paper.selected_gamma`: selected gamma.
- `paper.raw_partition`: raw consensus partition at the selected gamma.
- `paper.partition`: exact partition used in the paper.
- `paper.excluded_nodes`: one-based node indices excluded from the paper partition.
- `paper.equivalent_indices`: gamma indices producing an equivalent retained partition.

The selected-partition CSVs provide a lightweight view with columns:

```text
node_id,raw_community,community,included
```
