%% Load and inspect the exact structural partitions used in the paper

rootdir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootdir,'matlab'));
setup_paths(rootdir);

cfg_comm = cfg.default_community_config(rootdir);
cfg_comm.out_dir = fullfile(rootdir,'results','community_detection');
cfg_comm.networks = {'caliber','MTsat','gratio','delay','rate'};
cfg_comm.mode = 'precomputed';
cfg_comm.selection_mode = 'paper';
cfg_comm.plot_diagnostics = true;

results = analysis.run_community_detection([],cfg_comm);

% Example access:
% caliber_partition = results.caliber.selected.partition;
% mtsat_gamma = results.MTsat.selected.selected_gamma;
