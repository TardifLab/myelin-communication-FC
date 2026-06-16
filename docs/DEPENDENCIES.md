# Dependencies

## MATLAB

The package requires MATLAB and the Statistics and Machine Learning Toolbox. Functions used from that toolbox include Spearman correlation, F-distribution probabilities, and categorical dummy-variable construction.

The code uses modern MATLAB language and graphics features such as string arrays, `tiledlayout`, `boxchart`, and `exportgraphics` when available. A recent MATLAB release is recommended.

## Brain Connectivity Toolbox

Add the Brain Connectivity Toolbox to the MATLAB path separately. The manuscript analysis was developed with the `2019_03_03_BCT` release.

Required communication/community functions include:

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

Example setup:

```matlab
addpath(genpath('/path/to/2019_03_03_BCT'));
```

Because `navigation_wu` returns several path-length outputs, use a version with the signature:

```matlab
[sr, PL_bin, PL_wei, PL_dis, paths] = navigation_wu(L, D, max_hops)
```

## Eigenvector matching

The matched-eigenvector spectral analysis requires either:

- MATLAB's `matchpairs` function; or
- a compatible `munkres` implementation on the MATLAB path.

## Optional parallel execution

The Parallel Computing Toolbox is optional. It is used only when:

```matlab
cfg_comm.use_parallel = true;
```

for community-detection reruns. Precomputed paper partitions and all other default workflows run serially.

## Included helpers

The repository includes:

- `fcn_randz.m` for z-Rand partition stability;
- communication-model wrappers;
- lightweight plotting compatibility helpers;
- the full `+bicsdual` regression suite.

No Python or R environment is required for the supported MATLAB workflows.

## Optional manuscript-era visualization helpers

The supported CSV/MAT workflows do not require lab-specific surface plotting. The advanced function `bicsdual.sa_corr_nodal` can optionally call `plot_conn_surf` and `setsurf`; those functions are not bundled and are not used by `examples/run_full_replication.m`.
