%% Run the paper spectral analysis on the bundled Schaefer-400 data

rootdir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootdir,'matlab'));
setup_paths(rootdir);

% Add the Brain Connectivity Toolbox to the MATLAB path before running.
cfg_data = cfg.default_config();
[mats,pinfo] = io.load_inputs(cfg_data);
cm = analysis.compute_communication_models(mats,cfg_data);

cfg_spec = cfg.default_spectral_config(rootdir);
cfg_spec.out_dir = fullfile(rootdir,'results','spectral_analysis');
cfg_spec.matched_indices_mode = 'paper';

results = analysis.run_spectral_analysis(mats,cm,pinfo,cfg_spec); %#ok<NASGU>
