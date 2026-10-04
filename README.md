# Myelin-sensitive communication and functional connectivity

MATLAB code and analysis-ready derivative data for analyses examining how microstructure-weighted connectomes drive patterns of network communication and relate to multimodal functional connectivity.

The supported public workflows cover:

1. structural and functional matrix loading and edge-table export;
2. six network communication models;
3. communication-model correlations across all, structurally connected, and structurally disconnected node pairs;
4. hierarchical and dual-BICS nested regression at global, RSN-pair, and node levels;
5. structural-network community detection, including exact precomputed paper partitions and an optional two-pass rerun workflow;
6. structural spectral alignment, communication spectral fingerprints, and matched-eigenvector topology analyses.

The original messy scripts are retained under `legacy_scripts/` for provenance. New users should use the configuration-driven functions and examples described below.

## Repository structure

```text
data/
  schaefer-400/                  analysis-ready 400-node derivative matrices
  community_detection/          precomputed gamma sweeps and paper partitions
  spectral/                     paper matched-eigenvector index selection

docs/                            detailed workflow and output documentation
examples/
  run_full_replication.m         complete supported replication workflow
  run_minimal_example.m          small synthetic smoke test
  run_community_detection_example.m
  rerun_community_detection_example.m
  run_spectral_analysis_example.m
  run_spectral_smoke_test.m

matlab/
  setup_paths.m
  run_all.m                       convenience wrapper for one configured regression run
  +cfg/                           default configurations
  +io/                            input loading and edge-table export
  +prep/                          network normalization and delay/rate preparation
  +analysis/                      user-facing workflow entry points
  +bicsdual/                      nested regression framework
  +community/                     community-detection workflow
  +spectral/                      spectral workflow and outputs
  communication/                 communication-model implementations
  topology_spectra/              spectral/topology computation helpers

legacy_scripts/                  original messy analysis scripts
results/                         generated outputs; ignored by Git
```

## Requirements

- MATLAB with the Statistics and Machine Learning Toolbox.
- Brain Connectivity Toolbox functions used by the analyses. The manuscript workflow was developed with the BCT `2019_03_03` release.
- For matched-eigenvector assignment, either MATLAB's `matchpairs` function or a compatible `munkres` implementation.
- Parallel Computing Toolbox is optional and used only when explicitly enabled for community-detection reruns.

Required BCT functions include:

```text
distance_wei_floyd
navigation_wu
search_information
path_transitivity
diffusion_efficiency
community_louvain
agreement
consensus_und
```

See [`docs/DEPENDENCIES.md`](docs/DEPENDENCIES.md) for details.

## Setup

Clone or download the repository, then start MATLAB and run:

```matlab
rootdir = '/path/to/myelin-communication-FC';

% Add the Brain Connectivity Toolbox separately.
addpath(genpath('/path/to/2019_03_03_BCT'));

% Make setup_paths visible, then add the repository code.
addpath(fullfile(rootdir, 'matlab'));
setup_paths(rootdir);
```

The default configuration points to the derivative data already bundled under `data/schaefer-400/`:

```matlab
cfg_data = cfg.default_config();
[mats, pinfo] = io.load_inputs(cfg_data);
```

No source-code editing is required to run the included data. For custom data, override `cfg_data.input.*` after calling `cfg.default_config()`.

## Recommended full replication workflow

Open:

```text
examples/run_full_replication.m
```

Edit only the user settings at the top, particularly `bct_dir` and `output_root`, then run the script:

```matlab
run(fullfile(rootdir, 'examples', 'run_full_replication.m'));
```

The script performs the following supported workflow.

### Step 1: Load and validate bundled derivatives

The package loads five structural inputs and eight FC targets using a common Schaefer-400 node ordering:

- tract caliber;
- MTsat myelin density;
- tract-specific g-ratio;
- tract length;
- Euclidean distance;
- in-sample BOLD FC;
- out-of-sample BOLD FC;
- MEG delta, theta, alpha, beta, low-gamma, and high-gamma FC.

`nodes.csv` provides the RSN label for each node. It is required for meaningful RSN-pair regression summaries. The sensory-association axis used by the spectral module is stored in `sa_axis.mat`.

### Step 2: Export an edge table

```matlab
edges = io.mk_edgescsv_from_mat(cfg_data, output_csv);
```

The table contains one row per unique undirected node pair, structural variables, all FC targets, and the RSN labels of both nodes.

### Step 3: Compute communication models

```matlab
cm = analysis.compute_communication_models(mats, cfg_data);
```

The communication-model order is:

1. `SPE` — shortest-path efficiency
2. `NE` — navigation efficiency
3. `SIE` — inverse search information
4. `PT` — path transitivity
5. `CMY` — communicability
6. `DE` — diffusion efficiency

Raw communication matrices are retained for descriptive correlations and spectral fingerprints. Regression predictors are log-transformed when required by the manuscript workflow and then z-scored. FC targets and Euclidean distance are also z-scored for regression.

The communication summary reports Spearman correlations separately for:

- `all` — all unique node pairs;
- `con` — pairs connected in the empirical structural connectome;
- `dis` — structurally disconnected pairs.

### Step 4: Run main nested regressions

The full replication script processes one route-diffusion pair per call and loops over the four paper pairs:

```text
SPE-CMY
SPE-DE
NE-CMY
NE-DE
```

The main analyses add MTsat-, g-ratio-, and delay-sensitive communication predictors as one combined block:

```matlab
cfg_run.myelin_predictors = {'MTsat','gratio','delay'};
cfg_run.model_modes = 'all';
cfg_run.model_levels = [1 2];
```

`model_levels = [1 2]` produces:

- `L1_route` — routing-side communication model;
- `L1_diff` — diffusion-side communication model;
- `L2_both` — both communication regimes together.

The five tested regression conditions are:

| Output condition | `interaction` | `interact_at` | `effect_mode` | Interpretation |
|---|---|---|---|---|
| `main-effect` | `none` | `both` | `incremental` | Added myelin-sensitive main effects |
| `interact-caliber` | `caliber` | `both` | `incremental` | Added myelin × caliber terms |
| `interact-ed` | `ed` | `both` | `incremental` | Added myelin × Euclidean-distance terms |
| `interact-both` | `both` | `both` | `incremental` | Added both interaction families |
| `total-contribution` | `both` | `both` | `total` | Total myelin-related contribution, including main and interaction terms |

`total` is an `effect_mode`, not an `interaction` value. This is the correct total-contribution configuration:

```matlab
cfg_run.interaction = 'both';
cfg_run.interact_at = 'both';
cfg_run.effect_mode = 'total';
```

Nested-regression significance is assessed with an analytic nested F-test. Upper-tail probabilities are computed directly; values below floating-point resolution are stored as MATLAB's `realmin` rather than literal zero.

#### Multiple-comparison correction for the main regression analyses

`examples/run_full_replication.m` now saves raw Stage A results for all four
communication-model pairs and all five conditions, even when Stage B is
disabled. It then applies Benjamini–Hochberg (BH) FDR using
`bics_fdr_saved_results`. The settings at the top of the script are:

```matlab
paper_myelin_predictors = {'MTsat','gratio','delay'};
run_fdr = true;
fdr_alpha = 0.05;
```

Remove `'delay'` for the two-predictor sensitivity analysis. The main-regression
directory includes the selected predictor names and a unique run identifier.
Each invocation starts a new main-regression run; it does not automatically
resume a partially completed run. Other replication output directories retain
their existing behavior.

Correction families are separate for each experimental condition, each spatial
scale, and individual versus combined communication models. Each family
includes every analyzed FC target and location, irrespective of the subset
displayed in the main figures. With the paper's eight FC targets, seven RSNs,
and 400 nodes, the families are:

| Spatial scale | Individual models | Routing–diffusion combinations |
|---|---|---|
| Global | 4 × 8 = 32 | 4 × 8 = 32 |
| RSN | 4 × 8 × 28 = 896 | 4 × 8 × 28 = 896 |
| Node | 4 × 8 × 400 = 12,800 | 4 × 8 × 400 = 12,800 |

Individual models are SPE, NE, CMY, and DE. Combined models are SPE–CMY,
SPE–DE, NE–CMY, and NE–DE. Repeated individual-model outputs across pair files
must agree exactly and are counted only once. RSN hypotheses include the seven
within-network and 21 between-network blocks; mirrored entries are not counted
twice. The main workflow uses `model_modes='all'`, `model_levels=[1 2]`, and
`use_lower_only=true`. The correction helper rejects multiple predictor-menu
specifications (`M>1`); the supplementary `model_modes='single'` analyses remain
raw and require a separate, explicitly defined correction plan.

The rejection rule is adjusted p < `fdr_alpha`. Missing p-values occupy a family
slot as p=1 internally, retain NaN adjusted values, and are not rejected. Stored
zero p-values are accepted with a warning; FDR cannot recover numerical
precision lost in an earlier fit. FDR also does not repair violations of the
assumptions underlying the supplied analytic F-test p-values.

The `FDR` directory contains:

- `BICS_FDR_results.mat`: structure `F`, retaining raw effects/p-values and
  adding `Q_all`/`H_all` and `global_q`, `ntwk_q`, `node_q` with rejection masks.
- `BICS_FDR_family_summary.csv`: counts and rejections for every family.
- `BICS_FDR_global_tests.csv`: global effects, raw p, adjusted p, and rejection.
- `BICS_FDR_global_<condition>.csv`, `BICS_FDR_network_<condition>.csv`, and
  `BICS_FDR_node_<condition>.csv`: one copy of each unique model/location test,
  with raw `deltaR2` and `p`, adjusted `q`, `reject_fdr`, and the display-only
  `deltaR2_fdr` column. Failed/missing tests are NaN in the display column.

Raw estimates are never overwritten. Obtain separate analysis/display copies:

```matlab
z = load(fullfile(main_out,'FDR','BICS_FDR_results.mat'),'F');
F = z.F;
[Delta_raw,stats_raw] = bics_fdr_view(F,1,'total','raw');
[Delta_display,stats_display] = bics_fdr_view(F,1,'total','fdr');
```

The pair order is SPE–CMY, SPE–DE, NE–CMY, NE–DE. Use the FDR view for
thresholded bars, RSN heatmaps and surfaces. Non-rejected values are masked with
NaN, not zero; plotting code must render missing values explicitly. Use the
raw view for descriptive block averages, S–A correlations, and spin tests.
The helper changes which estimates are displayed, not their retained values.

To repeat correction without rerunning the fits:

```matlab
F = bics_fdr_saved_results(main_out, ...
    'FCLabels',mats.fc_labels,'RSNLabels',pinfo.clabels_short, ...
    'OutputDir',fullfile(main_out,'FDR_recheck'),'Overwrite',false);
analysis.write_fdr_tables(F,fullfile(main_out,'FDR_recheck'));
```

Existing correction outputs are protected. Choose a new output directory for
each recheck. Run `test_bics_fdr_saved_results` after `setup_paths` for synthetic
checks of BH answers, family sizes, duplicate handling, raw preservation,
display masks, unique-model CSV export, and invalid-input guards.

**Spin inference is separate.** The nested-regression `node_q` fields are not
adjusted spin p-values. `sa_corr_nodal` uses ordinary correlation p-values and
does not implement spatial null inference. Run the manuscript's separate
hemisphere-constrained spin workflow on unfiltered nodal effects and correct
those spin p-values using its own declared families before adding significance
markers. This FDR integration does not add or rerun spin tests.

### Step 5: Run supplemental individual-predictor models

The replication script separately runs:

```matlab
cfg_run.model_modes = 'single';
```

This tests MTsat, g-ratio, and delay predictors individually. Use only one of `single`, `pairs`, or `all` per run so model labels and outputs remain unambiguous.

### Step 6: Load the paper community partitions

The default reproducibility workflow uses the bundled precomputed gamma sweeps and exact final structural partitions:

```matlab
cfg_comm = cfg.default_community_config(rootdir);
cfg_comm.mode = 'precomputed';
cfg_comm.selection_mode = 'paper';
results = analysis.run_community_detection([], cfg_comm);
```

This is deterministic and does not rerun Louvain optimization. An optional two-pass rerun workflow is available for exploration or custom matrices:

```matlab
cfg_comm.mode = 'rerun';
cfg_comm.selection_mode = 'none';
results = analysis.run_community_detection(mats, cfg_comm);
```

The rerun automatically performs the coarse sweep, selects fine-sweep bounds, computes repeated partitions and consensus solutions, and exports stability diagnostics. Final partition selection remains explicit because high mean and low variance of z-Rand stability do not always identify the cleanest partition automatically.

See [`docs/COMMUNITY_DETECTION.md`](docs/COMMUNITY_DETECTION.md).

### Step 7: Run the spectral analysis

```matlab
cfg_spec = cfg.default_spectral_config(rootdir);
results = analysis.run_spectral_analysis(mats, cm, pinfo, cfg_spec);
```

The comprehensive spectral module includes:

- direct structural spectral alignment to caliber;
- communication-model spectral fingerprints;
- globally banded matched-eigenvector comparisons;
- selected-eigenvector topology summaries.

The caliber adjacency and normalized-Laplacian bases are computed directly. The paper-selected matched eigenvectors are loaded from:

```text
data/spectral/paper_matched_eigenvector_indices.csv
```

The Schaefer-400 sensory-association axis is loaded automatically from:

```text
data/schaefer-400/sa_axis.mat
```

Lightweight built-in MATLAB plots are shown by default and can optionally be saved.

See [`docs/SPECTRAL_ANALYSIS.md`](docs/SPECTRAL_ANALYSIS.md).

## Output organization

The full replication script writes to:

```text
results/paper_replication/
  edges_schaefer400.csv
  communication/
  regression/
    main_results_MTsat_gratio_delay/
      run_<timestamp>_<unique-id>/
        results_total_tmpSPE-CMY.mat
        results_total_tmpSPE-DE.mat
        results_total_tmpNE-CMY.mat
        results_total_tmpNE-DE.mat
        global/<condition>/
        network/<condition>/
        node/<condition>/
        FDR/
    supplementary_results/
      global/main-effect/
      network/main-effect/
      node/main-effect/
  community_detection/
  spectral_analysis/
```

Regression filenames encode the communication pair, myelin predictor grouping, and included model levels, for example:

```text
bics_deltaR2_SPE_DE_all_levels_L1L2.csv
```

The directory hierarchy encodes spatial scale and regression condition. Each CSV also contains route-model, diffusion-model, communication-pair, and condition metadata columns.

See [`docs/OUTPUTS.md`](docs/OUTPUTS.md).

## Single configured regression run

`run_all(cfg)` remains available as a convenience wrapper for one communication pair, predictor mode, and regression condition:

```matlab
cfg_run = cfg.default_config();
cfg_run.out_dir = fullfile(rootdir, 'results', 'single_run');
cfg_run.dual_bics_pairs = [1 6];
cfg_run.model_modes = 'all';
cfg_run.interaction = 'none';
cfg_run.effect_mode = 'incremental';
results = run_all(cfg_run);
```

For replication of all paper conditions and pairs, use `examples/run_full_replication.m` instead.

`run_all` returns raw regression results: one pair alone cannot supply the
paper's across-model FDR families. Use the full replication workflow to apply
that correction. The raw RSN CSV writer exports only one triangle including
within-network blocks for the symmetric undirected analysis.

## Custom data

All required connectivity matrices must be square, numeric, finite, symmetric, and use the same node ordering. Structural edge-weight matrices should share the same empirical binary support. Euclidean distance and FC matrices may be dense.

See [`docs/DATA_FORMAT.md`](docs/DATA_FORMAT.md).

## Legacy scripts and plotting

The supported workflows emphasize reusable CSV and MAT outputs plus lightweight diagnostic plots. Manuscript-specific plotting code, local path assumptions, and exploratory alternatives remain in `legacy_scripts/` for provenance and are not the recommended entry points.

## Citation and license

Please see `CITATION.cff` for manuscript citation info. This repository is distributed under the GNU General Public License v3.0; see `LICENSE`.
