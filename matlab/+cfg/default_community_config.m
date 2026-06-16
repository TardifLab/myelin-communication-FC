function cfg = default_community_config(rootdir)
%DEFAULT_COMMUNITY_CONFIG Configuration for structural-network community analysis.
%
% The default mode loads the precomputed gamma sweeps and exact partitions
% used in the paper. Set cfg.mode = 'rerun' to repeat the two-pass Louvain
% analysis on supplied or custom matrices.

if nargin < 1 || isempty(rootdir)
    rootdir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
end

cfg = struct();
cfg.root_dir = rootdir;
cfg.out_dir = fullfile(rootdir, 'results', 'community_detection');

% Networks to process. Accepted aliases include myelin/mtsat and g-ratio/gratio.
cfg.networks = {'caliber','MTsat','gratio','delay','rate'};

% 'precomputed': load the bundled second-pass gamma sweeps.
% 'rerun': repeat coarse and fine Louvain sweeps on input matrices.
cfg.mode = 'precomputed';
cfg.precomputed_dir = fullfile(rootdir, 'data', 'community_detection', 'precomputed');

% Partition selection after the fine sweep:
% 'paper'  = exact bundled paper partition (precomputed mode only)
% 'none'   = retain all candidates without selecting one
% 'score'  = maximize alpha*mean_zRand + (1-alpha)*(1-var_zRand)
% 'gamma'  = choose the sampled gamma nearest cfg.selected_gamma
cfg.selection_mode = 'paper';
cfg.selected_gamma = [];
cfg.score_alpha = 0.30;
cfg.min_community_size = 1;   % applied only to score/gamma selections

% Two-pass sweep settings used by rerun mode.
cfg.coarse.gamma_bounds = [0.1 5];
cfg.coarse.n_gamma = 100;
cfg.coarse.lower_count_threshold = 2;
cfg.coarse.upper_count_fraction = 0.5;
cfg.fine.n_gamma = 200;
cfg.fine.n_repetitions = 100;
cfg.fine.padding = 2;
cfg.fine.consensus_repetitions = 10;

% Reproducibility and execution.
cfg.random_seed = 13;
cfg.use_parallel = false;
cfg.symmetry_mode = 'average';  % 'average' or 'upper'

% Output and display.
cfg.export_csv = true;
cfg.save_results_mat = true;
cfg.plot_diagnostics = true;
cfg.save_figures = false;

% Optional direct matrix overrides for rerun mode, e.g.:
% cfg.network_matrices.MTsat = your_400_by_400_matrix;
cfg.network_matrices = struct();

% Delay/rate constants used when these matrices are derived from length and
% g-ratio rather than supplied directly.
cfg.delay.k = 6.0;
cfg.delay.axon_diameter_um = 2.5;
end
