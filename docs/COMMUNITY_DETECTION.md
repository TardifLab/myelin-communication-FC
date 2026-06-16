# Community detection

## Two supported modes

### Exact paper partitions

The default and recommended replication mode loads the standardized second-pass sweeps and exact final partitions used in the paper:

```matlab
cfg_comm = cfg.default_community_config(rootdir);
cfg_comm.mode = 'precomputed';
cfg_comm.selection_mode = 'paper';
results = analysis.run_community_detection([], cfg_comm);
```

Available structural networks are:

```text
caliber
MTsat
gratio
delay
rate
```

This mode is deterministic and does not rerun Louvain optimization.

### Rerun mode

```matlab
[mats,~] = io.load_inputs(cfg.default_config());

cfg_comm = cfg.default_community_config(rootdir);
cfg_comm.mode = 'rerun';
cfg_comm.networks = {'MTsat'};
cfg_comm.selection_mode = 'none';
results = analysis.run_community_detection(mats, cfg_comm);
```

Rerun mode:

1. transforms the selected structural network to positive Louvain weights;
2. runs a coarse gamma sweep;
3. derives the fine-sweep bounds automatically from community counts;
4. runs repeated Louvain optimizations at each fine gamma;
5. computes agreement and consensus partitions;
6. calculates mean and variance of pairwise z-Rand stability;
7. exports all candidate partitions and diagnostic metrics.

## Final selection

Available selection modes are:

```text
paper   exact bundled paper partition; precomputed mode only
none    do not select a final partition
score   automated mean/variance stability score
gamma   sampled gamma nearest cfg.selected_gamma
```

The default rerun mode is `none`. Automated stability-score selection is retained as a candidate method, but final interpretation may require visual review because the maximum score did not always produce the cleanest community structure in the manuscript analysis.

## Reproducibility

Community detection is stochastic. Rerun mode controls random seeds, but exact partitions can still depend on software environment and dependency versions. Use precomputed mode to reproduce the reported partitions exactly.

## Outputs

Each processed network receives its own directory containing:

```text
gamma_summary.csv
selected_partition.csv       when a partition is selected
community_results.mat
stability_diagnostics.png    only when save_figures=true
```

The bundled precomputed schema is documented in [`COMMUNITY_DATA_SCHEMA.md`](COMMUNITY_DATA_SCHEMA.md).
