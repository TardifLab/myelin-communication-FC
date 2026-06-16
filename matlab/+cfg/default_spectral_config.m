function config = default_spectral_config(rootdir)
%DEFAULT_SPECTRAL_CONFIG Configuration for the spectral workflow.
%
% The workflow combines:
%   1) Direct caliber-versus-target spectral comparison (legacy Module 1)
%   2) Communication-model spectral fingerprints (legacy Module 3)
%   3) Global-banded matched-eigenvector topology comparisons (Module 3)
%
% Existing functions under matlab/topology_spectra are reused for the
% fingerprint and matched-eigenvector computations.

if nargin < 1 || isempty(rootdir)
    rootdir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

config = struct();
config.root_dir = rootdir;
config.out_dir = fullfile(rootdir, 'results', 'spectral_analysis');

% Components to run.
config.run_alignment = true;
config.run_fingerprints = true;
config.run_matched_eigenvectors = true;

% Structural datasets. Caliber is the reference basis. The dataset called
% "delay" below is the delay-derived communication-rate structural matrix,
% matching the naming used in the manuscript analysis code.
config.reference_dataset = 'caliber';
config.target_datasets = {'MTsat','gratio','delay'};
config.dataset_labels = {'caliber','MTsat','gratio','delay'};
config.target_labels = {'MTsat','g-ratio','delay'};

% Module 1: direct pairwise comparison.
config.alignment_kmax = 50;
config.sign_align = true;

% Module 3: spectral basis and fingerprint settings.
config.kmax = 399;
config.windows = 'auto';
config.quantile_edges = [0 .05 .25 .55 1];
config.start_at = 2;
config.reference_dataset_index = 1;
config.communication_models = {'SPE','NE','SIE','PT','CMY','DE'};
config.operator_for_model = {'adj','adj','adj','adj','adj','lap'};
config.de_use_dhalf = true;
config.renormalize_covered = true;
config.deoverlap_windows = false;
config.split_dc_figure = true;
config.show_both_dc = true;

% Matched-eigenvector analysis.
config.match_mode = 'global-banded';
config.match_band = 3;
config.match_k_per_page = 399;
config.drop_laplacian_dc = true;
config.sa_axis = [];              % Optional explicit N x 1 sensory-association axis
config.sa_axis_file = fullfile(rootdir, 'data', 'schaefer-400', ...
    'sa_axis.mat');                 % Loaded automatically when sa_axis is empty
config.sa_corr_type = 'pearson';
config.compute_moran = true;

% Delay/rate reconstruction fallback.
config.delay = struct('k',6.0,'axon_diameter_um',2.5);

% Eigenvector subset used in the supplementary comparison.
% Modes: 'paper', 'custom', or 'all'.
config.matched_indices_mode = 'paper';
config.matched_indices_csv = fullfile(rootdir, 'data', 'spectral', ...
    'paper_matched_eigenvector_indices.csv');
config.custom_indices.adj = [];
config.custom_indices.lap = [];

% Lightweight plotting.
config.make_plots = true;
config.figure_visible = 'on';     % 'on' or 'off'
config.plot_selected_metrics = true;
% Selected-eigenvector plot: 'windowed', 'pooled', or 'both'.
% 'windowed' reproduces the supplemental-style per-window summary.
config.selected_metrics_plot_mode = 'windowed';
% Windowed bars: 'norm' reproduces the legacy relative scaling; 'raw'
% retains the original target-minus-caliber metric units.
config.selected_metrics_bar_mode = 'norm';
config.selected_metric_names = { ...
    'delta_SAr','delta_smooth','delta_PR','delta_modR2','delta_MoranI'};

% The reused robust matching helper creates large diagnostic figures.
% Keep these hidden by default; the add-on creates compact plots instead.
config.show_matched_helper_plots = false;

% Optional legacy-style all-k boxcharts. This can create many figures.
config.plot_windowed_topology_boxcharts = false;

config.save_figures = false;
config.figure_format = 'png';

% Output controls.
config.write_csv = true;
config.save_mat = true;
config.verbose = true;
end
