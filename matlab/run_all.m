function results = run_all(config)
%RUN_ALL Run one configured communication and nested-regression analysis.
%
% This convenience wrapper loads inputs, computes communication models,
% runs one configured route-diffusion pair and regression condition, and
% writes analysis-ready outputs. Use examples/run_full_replication.m to loop
% over all paper model pairs and conditions.
%
% Usage:
%   setup_paths
%   config = cfg.default_config();   % uses bundled Schaefer-400 data
%   results = run_all(config);

if nargin < 1 || isempty(config)
    config = cfg.default_config();
end
setup_paths(fileparts(fileparts(mfilename('fullpath'))));
if ~exist(config.out_dir,'dir'), mkdir(config.out_dir); end

[mats, pinfo] = io.load_inputs(config);
cm = analysis.compute_communication_models(mats, config);
analysis.write_communication_tables(cm, mats, config);

results = analysis.run_nested_regression(cm, mats, pinfo, config);
save(fullfile(config.out_dir,'paper_results.mat'), 'results', 'cm', 'mats', 'pinfo', '-v7.3');
analysis.write_bics_tables(results, config);

fprintf('Done. Results written to: %s\n', config.out_dir);
end
