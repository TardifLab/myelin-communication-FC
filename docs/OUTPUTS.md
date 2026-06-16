# Output files

## Full replication

`examples/run_full_replication.m` writes to a user-selected output root, by default:

```text
results/paper_replication/
```

## Edge and communication outputs

```text
edges_schaefer400.csv
communication/
  communication_edge_correlations.csv
  communication_models.mat
```

The edge correlation table includes `all`, `con`, and `dis` edge groups.

## Regression outputs

```text
regression/
  main_results/
    global/<condition>/
    network/<condition>/
    node/<condition>/
  supplementary_results/
    global/main-effect/
    network/main-effect/
    node/main-effect/
```

Conditions are:

```text
main-effect
interact-caliber
interact-ed
interact-both
total-contribution
```

Filenames encode communication pair, predictor mode, and levels:

```text
bics_deltaR2_<route>_<diffusion>_<mode>_levels_L1L2.csv
```

Examples:

```text
bics_deltaR2_SPE_CMY_all_levels_L1L2.csv
bics_deltaR2_NE_DE_single_levels_L1L2.csv
```

Global CSV columns include:

```text
route_model,diffusion_model,communication_pair,condition,
level,model,fc,deltaR2,p
```

Network and node tables add RSN-pair or node identifiers.

## Community outputs

```text
community_detection/
  paper_partitions/<network>/
    gamma_summary.csv
    selected_partition.csv
    community_results.mat
```

Optional reruns are written separately.

## Spectral outputs

```text
spectral_analysis/
  alignment/
  fingerprints/
    spectral_fingerprints.csv
    spectral_windows.csv
  matched_eigenvectors/
    matched_eigenvector_metrics_all.csv
    matched_eigenvector_metrics_selected.csv
    selected_eigenvector_indices_used.csv
  reference_basis_eigenvalues.csv
  spectral_results.mat
  figures/                         when save_figures=true
```

## Generated files and Git

The `results/` directory is ignored by Git except for `.gitkeep`. Generated outputs should generally not be committed unless a specific release requires archived reference results.
