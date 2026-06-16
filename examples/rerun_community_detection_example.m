%% Rerun the two-pass community analysis on one structural network
%
% This is an exploratory stochastic rerun. To load the exact partitions
% used in the paper, use run_community_detection_example.m instead.

rootdir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootdir,'matlab'));
setup_paths(rootdir);

% Add the Brain Connectivity Toolbox to the MATLAB path before running.
main_cfg = cfg.default_config();
[mats,~] = io.load_inputs(main_cfg);

cfg_comm = cfg.default_community_config(rootdir);
cfg_comm.out_dir = fullfile(rootdir,'results','community_detection_rerun');
cfg_comm.networks = {'MTsat'};
cfg_comm.mode = 'rerun';
cfg_comm.selection_mode = 'none'; % review diagnostics before choosing gamma
cfg_comm.use_parallel = false;

results = analysis.run_community_detection(mats,cfg_comm); %#ok<NASGU>

% After reviewing the stability diagnostics, select a sampled gamma:
% cfg_comm.selection_mode = 'gamma';
% cfg_comm.selected_gamma = 2.04;
% results = analysis.run_community_detection(mats,cfg_comm);
