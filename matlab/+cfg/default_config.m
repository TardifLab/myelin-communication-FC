function cfg = default_config()
%DEFAULT_CONFIG Default configuration for the bundled paper derivatives.
%
% The repository includes analysis-ready Schaefer-400 structural and
% functional connectivity matrices under data/schaefer-400. This function
% points to those files by default so the public workflow can run without
% editing source code. For custom data, override cfg.input.* after calling
% this function.
%
% Nested-regression predictor grouping:
%   'single' = test each listed myelin predictor separately
%   'pairs'  = test pairwise combinations
%   'all'    = test all listed predictors as one combined block
%
% Use one model_modes value per run. The main-paper analyses use 'all'; the
% individual-predictor supplemental analyses use 'single'.

rootdir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir = fullfile(rootdir, 'data', 'schaefer-400');

cfg = struct();
cfg.root_dir = rootdir;
cfg.out_dir  = fullfile(rootdir, 'results');

% -------------------------------------------------------------------------
% Bundled input matrices. All matrices are 400 x 400, symmetric, and share
% the same Schaefer-400 node ordering.
% -------------------------------------------------------------------------
cfg.input = struct();
cfg.input.caliber_mat   = fullfile(data_dir, 'caliber.mat');
cfg.input.mtsat_mat     = fullfile(data_dir, 'myelin_density.mat');
cfg.input.gratio_mat    = fullfile(data_dir, 'gratio_tract_specific.mat');
cfg.input.length_mat    = fullfile(data_dir, 'length.mat');
cfg.input.euclidean_mat = fullfile(data_dir, 'euclidean_distance.mat');

cfg.input.fc_mats = {
    fullfile(data_dir, 'fc_boldin.mat')
    fullfile(data_dir, 'fc_boldout.mat')
    fullfile(data_dir, 'fc_meg_delta.mat')
    fullfile(data_dir, 'fc_meg_theta.mat')
    fullfile(data_dir, 'fc_meg_alpha.mat')
    fullfile(data_dir, 'fc_meg_beta.mat')
    fullfile(data_dir, 'fc_meg_gammalo.mat')
    fullfile(data_dir, 'fc_meg_gammahi.mat')
};

cfg.input.fc_labels = {
    'BOLDin'
    'BOLDout'
    'delta'
    'theta'
    'alpha'
    'beta'
    'gammalo'
    'gammahi'
};

% Node metadata are required for meaningful RSN-pair summaries.
cfg.input.nodes_csv = fullfile(data_dir, 'nodes.csv');

% -------------------------------------------------------------------------
% Communication and nested-regression options
% -------------------------------------------------------------------------
cfg.parcellation = 'schaefer-400';
cfg.length_transform = 'log';
cfg.communication_models = {'SPE','NE','SIE','PT','CMY','DE'};

% One route-diffusion pair per regression run. Examples:
%   [1 5] SPE-CMY, [1 6] SPE-DE, [2 5] NE-CMY, [2 6] NE-DE
cfg.dual_bics_pairs = [1 5];

cfg.myelin_predictors = {'MTsat','gratio','delay'};
cfg.model_modes = 'all';
cfg.model_levels = [1 2];

% interaction: 'none', 'caliber', 'ed', or 'both'
% interact_at: 'route', 'diff', or 'both'
% effect_mode: 'incremental' or 'total'
cfg.interaction = 'none';
cfg.interact_at = 'both';
cfg.effect_mode = 'incremental';

cfg.use_lower_only = true;
cfg.alpha = 0.05;

% Delay proxy constants used in the manuscript analysis.
cfg.delay.k = 6.0;                 % m/s per micron
cfg.delay.axon_diameter_um = 2.5;

% Stage B permutation/sensitivity analysis is optional and substantially
% slower than the Stage A nested regressions.
cfg.run_stageB = false;
cfg.stageB.nperm = 200;
cfg.stageB.eps_scale = 0.05;
cfg.stageB.rng_seed = 13;
cfg.stageB.perm_mode = 'yresid';
cfg.stageB.elasticity_mode = 'zrow';
end
