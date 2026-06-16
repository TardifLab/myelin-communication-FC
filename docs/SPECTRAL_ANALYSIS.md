# Spectral analysis

The comprehensive spectral workflow includes:

1. direct alignment of MTsat, g-ratio, and delay-derived rate structural spectra to caliber;
2. spectral fingerprints of six communication models in the caliber basis;
3. globally banded matching of caliber and myelin-sensitive eigenvectors;
4. topology comparisons for all matched eigenvectors and for the paper-selected subset.

The manually windowed intermediate workflow from the legacy script is not required. The caliber adjacency and normalized-Laplacian bases are computed directly.

## Basic use

```matlab
cfg_data = cfg.default_config();
[mats,pinfo] = io.load_inputs(cfg_data);
cm = analysis.compute_communication_models(mats,cfg_data);

cfg_spec = cfg.default_spectral_config(rootdir);
results = analysis.run_spectral_analysis(mats,cm,pinfo,cfg_spec);
```

The fingerprint analysis intentionally uses `cm.raw`, matching the manuscript analysis. Transformed/z-scored regression predictors are not used for fingerprints.

The structural dataset labeled `delay` in the spectral module is the delay-derived signaling-rate network used in the manuscript spectral comparisons.

## Existing helpers

The workflow reuses functions under `matlab/topology_spectra/`, including:

```text
normalize_operators.m
spectral_fingerprint_windows.m
spectral_alignment_windows_robust.m
rowStandardizeWeights.m
wins_for.m
```

The optional all-k topology summaries also use `delta_topology_boxcharts.m`.

## Paper matched-eigenvector subset

The default setting loads the indices used in the paper:

```matlab
cfg_spec.matched_indices_mode = 'paper';
```

The selection is stored in:

```text
data/spectral/paper_matched_eigenvector_indices.csv
```

Custom indices can be supplied with:

```matlab
cfg_spec.matched_indices_mode = 'custom';
cfg_spec.custom_indices.adj = [7 8 9 12];
cfg_spec.custom_indices.lap = [7 8 9 12];
```

To summarize all matched eigenvectors:

```matlab
cfg_spec.matched_indices_mode = 'all';
```

## Sensory-association axis

By default, the module automatically loads:

```text
data/schaefer-400/sa_axis.mat
```

The loader accepts a numeric vector of the expected node count. An explicit vector takes precedence:

```matlab
cfg_spec.sa_axis = custom_sa_axis(:);
```

or a different file can be supplied:

```matlab
cfg_spec.sa_axis_file = '/path/to/sa_axis.mat';
```

## Plot controls

Lightweight pairwise alignment, fingerprint, and selected matched-eigenvector plots are shown by default.

```matlab
cfg_spec.make_plots = true;
cfg_spec.figure_visible = 'on';
cfg_spec.save_figures = false;
```

The selected matched-eigenvector summary is windowed by default, matching the supplemental analysis: each panel is a topology metric, x-axis groups are spectral windows, and bars are target datasets.

```matlab
cfg_spec.selected_metrics_plot_mode = 'windowed'; % 'windowed', 'pooled', or 'both'
cfg_spec.selected_metrics_bar_mode = 'norm';      % 'norm' or 'raw'
```

`'norm'` scales each metric by its maximum absolute window mean within ADJ or LAP, reproducing the legacy relative-scale visualization. Use `'raw'` to retain the original target-minus-caliber metric units.

The reused global-banded matching helper contains larger diagnostic figures. These remain hidden unless:

```matlab
cfg_spec.show_matched_helper_plots = true;
```

Optional all-k topology boxcharts can be enabled with:

```matlab
cfg_spec.plot_windowed_topology_boxcharts = true;
```

## Outputs

```text
alignment/
fingerprints/
matched_eigenvectors/
reference_basis_eigenvalues.csv
spectral_results.mat
```

The matched-eigenvector tables contain reference and matched target indices, cosine similarity, topology metrics, and target-minus-caliber differences.
