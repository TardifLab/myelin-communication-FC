%% Minimal synthetic smoke test
%
% Generates a small synthetic dataset and runs one configured communication
% and nested-regression analysis. The Brain Connectivity Toolbox must be on
% the MATLAB path.

rootdir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(rootdir,'matlab'));
setup_paths(rootdir);

example_dir = fullfile(rootdir,'data','example');
if ~exist(example_dir,'dir'), mkdir(example_dir); end

rng(1);
N = 40;
A = rand(N); A = triu(A,1); A = A + A'; A(A < 0.75) = 0;
Dtg = A; save(fullfile(example_dir,'caliber.mat'),'Dtg');
Dtg = A .* (0.5 + rand(N)); Dtg = (Dtg+Dtg')/2; save(fullfile(example_dir,'mtsat.mat'),'Dtg');
Dtg = A .* (0.6 + 0.2*rand(N)); Dtg = (Dtg+Dtg')/2; save(fullfile(example_dir,'gratio.mat'),'Dtg');
coords = rand(N,3);
Dtg = squareform(pdist(coords)) .* 100; Dtg(A==0)=0; save(fullfile(example_dir,'edge_length.mat'),'Dtg');
Dtg = squareform(pdist(coords)); save(fullfile(example_dir,'euclidean_distance.mat'),'Dtg');
Dtg = corr(randn(200,N)); Dtg(1:N+1:end)=0; save(fullfile(example_dir,'fc_BOLDin.mat'),'Dtg');

node_id = (1:N)'; rsn_id = ceil((1:N)'/10); rsn = "RSN" + string(rsn_id);
writetable(table(node_id,rsn_id,rsn), fullfile(example_dir,'nodes.csv'));

config = cfg.default_config();
config.input.caliber_mat = fullfile(example_dir,'caliber.mat');
config.input.mtsat_mat = fullfile(example_dir,'mtsat.mat');
config.input.gratio_mat = fullfile(example_dir,'gratio.mat');
config.input.length_mat = fullfile(example_dir,'edge_length.mat');
config.input.euclidean_mat = fullfile(example_dir,'euclidean_distance.mat');
config.input.fc_mats = {fullfile(example_dir,'fc_BOLDin.mat')};
config.input.fc_labels = {'BOLDin'};
config.input.nodes_csv = fullfile(example_dir,'nodes.csv');
config.out_dir = fullfile(rootdir,'results','example');
config.dual_bics_pairs = [1 5];
config.model_modes = 'all';
config.model_levels = [1 2];
config.interaction = 'none';
config.effect_mode = 'incremental';

results = run_all(config); %#ok<NASGU>
