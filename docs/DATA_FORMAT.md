# Data format

## Bundled derivative data

The default configuration uses analysis-ready Schaefer-400 derivatives in:

```text
data/schaefer-400/
```

Files include:

| File | Contents |
|---|---|
| `caliber.mat` | tract-caliber structural edge weights |
| `myelin_density.mat` | MTsat myelin-density edge weights |
| `gratio_tract_specific.mat` | tract-specific g-ratio edge weights |
| `length.mat` | tract length on empirical structural edges |
| `euclidean_distance.mat` | dense inter-node Euclidean distance |
| `fc_boldin.mat` | in-sample BOLD FC |
| `fc_boldout.mat` | out-of-sample BOLD FC |
| `fc_meg_*.mat` | MEG FC by frequency band |
| `nodes.csv` | node IDs and resting-state-network labels |
| `sa_axis.mat` | Schaefer-400 sensory-association axis |

The matrices are 400 × 400 and share the same node ordering. The structural edge-weight matrices and tract length share the same binary support.

## Matrix files

`io.load_matrix_var` supports:

- `.mat` files containing `Dtg`;
- `.mat` files containing one square numeric matrix;
- square numeric `.csv` or `.txt` files.

All connectivity matrices used together must be:

- square and numeric;
- finite;
- symmetric, apart from negligible floating-point error;
- aligned to the same node order.

The diagonal should normally be zero.

## Node metadata

`nodes.csv` must contain one row per node. The minimal supported format is:

```text
node_id,rsn
1,Visual
2,Visual
...
```

An optional numeric `rsn_id` column may also be supplied:

```text
node_id,rsn_id,rsn
```

If node metadata are absent, global and node-level regression can still run, but all nodes are treated as one network and RSN-pair summaries are not meaningful.

## Functional targets

Provide one matrix per FC modality or frequency band and list them in the same order as their labels:

```matlab
cfg.input.fc_mats = {
    '/path/to/fc_bold.mat'
    '/path/to/fc_alpha.mat'
};

cfg.input.fc_labels = {
    'BOLD'
    'alpha'
};
```

Separate `.mat` files are recommended for MATLAB users rather than requiring an external NumPy loader.

## Custom configuration

Start from the defaults and override paths without editing package source files:

```matlab
cfg_custom = cfg.default_config();
cfg_custom.input.caliber_mat = '/path/to/caliber.mat';
cfg_custom.input.mtsat_mat = '/path/to/mtsat.mat';
cfg_custom.input.gratio_mat = '/path/to/gratio.mat';
cfg_custom.input.length_mat = '/path/to/length.mat';
cfg_custom.input.euclidean_mat = '/path/to/euclidean.mat';
cfg_custom.input.nodes_csv = '/path/to/nodes.csv';
```

For custom parcellations, also provide a matching sensory-association vector if the spectral `delta_SAr` metric is desired:

```matlab
cfg_spec.sa_axis = custom_sa_axis(:);
```
