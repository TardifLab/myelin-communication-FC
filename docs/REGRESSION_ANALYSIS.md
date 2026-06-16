# Nested regression

## Overview

The Stage A workflow tests whether myelin-sensitive communication predictors explain functional connectivity beyond lower-order structural and topological predictors. Analyses are performed separately at:

- global edge level;
- resting-state-network pair level;
- node level.

## Predictor preparation

`analysis.compute_communication_models` produces both raw and regression-ready communication matrices:

- `cm.raw`: untransformed communication estimates used for descriptive correlations and spectral fingerprints;
- `cm.z`: log-transformed when required by skewness and then z-scored;
- `cm.ALLX`: the transformed/z-scored predictors supplied to regression.

FC targets and normalized Euclidean distance are also z-scored before modeling.

## Communication-model pairs

The model order is:

```text
1 SPE
2 NE
3 SIE
4 PT
5 CMY
6 DE
```

The paper replication evaluates:

```matlab
model_pairs = [1 5; 1 6; 2 5; 2 6];
```

Each regression call should contain one route-diffusion pair. `write_bics_tables` encodes the selected pair in the output filename and metadata.

## Predictor grouping

Use one grouping mode per run:

```matlab
cfg.model_modes = 'single'; % individual predictors
cfg.model_modes = 'pairs';  % pairwise predictor blocks
cfg.model_modes = 'all';    % all listed predictors as one block
```

The main-paper workflow uses:

```matlab
cfg.myelin_predictors = {'MTsat','gratio','delay'};
cfg.model_modes = 'all';
```

The supplemental individual-predictor workflow uses `single` in a separate run.

## Model levels

```matlab
cfg.model_levels = [1 2];
```

produces:

- `L1_route`;
- `L1_diff`;
- `L2_both`.

By design, `L1_route` is identical between model pairs sharing the same routing model, and `L1_diff` is identical between pairs sharing the same diffusion model.

## Conditions

| Condition | Settings |
|---|---|
| Main effects | `interaction='none'`, `effect_mode='incremental'` |
| Myelin × caliber | `interaction='caliber'`, `interact_at='both'`, `effect_mode='incremental'` |
| Myelin × ED | `interaction='ed'`, `interact_at='both'`, `effect_mode='incremental'` |
| Both interaction families | `interaction='both'`, `interact_at='both'`, `effect_mode='incremental'` |
| Total contribution | `interaction='both'`, `interact_at='both'`, `effect_mode='total'` |

Do not set `interaction='total'`. `total` is an effect mode.

## Significance testing

Stage A uses an analytic nested F-test based on the numerical ranks of reduced and full models. Rank-deficient designs are solved with minimum-norm least squares. A constant or redundant predictor contributes no additional rank.

The F upper-tail probability is computed directly. If the probability underflows to zero in double precision, the saved value is floored at MATLAB's `realmin`.

Stage B permutation sensitivity and elasticity are separate optional analyses:

```matlab
cfg.run_stageB = true;
```

They are slower and are disabled in the default workflow.

## Output layout

```text
<regression output>/
  global/<condition>/
  network/<condition>/
  node/<condition>/
```

A typical filename is:

```text
bics_deltaR2_SPE_DE_all_levels_L1L2.csv
```

Each file includes `deltaR2`, `p`, FC labels, model labels, communication-pair metadata, and scale-specific identifiers.
