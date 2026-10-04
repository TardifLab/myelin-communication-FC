%% Paper full replication workflow
%
% This script runs the supported public workflows on the bundled
% Schaefer-400 derivative data:
%   1. Load and export analysis-ready edge data
%   2. Compute communication models and edgewise correlations
%   3. Run Stage A nested regressions for all tested model pairs/conditions
%      Save raw outputs and apply the paper's across-model BH-FDR families
%   4. Run supplemental individual-myelin-predictor regressions
%   5. Load the exact paper community partitions
%   6. Run the comprehensive spectral analysis
%
% Community-detection reruns and Stage B permutation/sensitivity analyses
% are available but disabled by default because they are computationally
% expensive. Edit the USER SETTINGS section as needed.

%% USER SETTINGS
rootdir = fileparts(fileparts(mfilename('fullpath')));

% Add the Brain Connectivity Toolbox separately before running. Either set
% this to its folder or leave it empty if BCT is already on the MATLAB path.
bct_dir = '';

output_root = fullfile(rootdir, 'results', 'paper_replication');

run_edge_export                    = true;
run_communication_models           = true;
run_main_regressions               = true;
run_supplementary_single_predictors = true;
run_precomputed_communities        = true;
run_community_rerun                = false;
run_spectral_analysis              = true;

% Optional expensive settings
run_stageB = false;
community_rerun_networks = {'MTsat'};

% Main-paper nested F-test correction. 'all' fits these predictors jointly.
paper_myelin_predictors = {'MTsat','gratio','delay'};
run_fdr = true;
fdr_alpha = 0.05;

%% SETUP
if ~isempty(bct_dir)
    if ~exist(bct_dir, 'dir')
        error('BCT directory not found: %s', bct_dir);
    end
    addpath(genpath(bct_dir));
end

addpath(fullfile(rootdir, 'matlab'));
setup_paths(rootdir);

required_bct = {};
if run_communication_models || run_main_regressions || ...
        run_supplementary_single_predictors || run_spectral_analysis
    required_bct = [required_bct; {
        'distance_wei_floyd'
        'navigation_wu'
        'search_information'
        'path_transitivity'
        'diffusion_efficiency'
    }];
end
if run_community_rerun
    required_bct = [required_bct; {
        'community_louvain'
        'agreement'
        'consensus_und'
    }];
end
required_bct = unique(required_bct, 'stable');
missing = required_bct(cellfun(@(f) exist(f,'file') ~= 2, required_bct));
if ~isempty(missing)
    error(['Missing required Brain Connectivity Toolbox functions: %s\n', ...
        'Add BCT to the MATLAB path and rerun this script.'], ...
        strjoin(missing, ', '));
end

if ~exist(output_root, 'dir'), mkdir(output_root); end

%% LOAD BUNDLED DATA
cfg_data = cfg.default_config();
cfg_data.out_dir = output_root;
[mats, pinfo] = io.load_inputs(cfg_data);

fprintf('Loaded %d nodes and %d FC datasets.\n', ...
    size(mats.caliber,1), numel(mats.FC));

%% EDGE TABLE
if run_edge_export
    edge_file = fullfile(output_root, 'edges_schaefer400.csv');
    io.mk_edgescsv_from_mat(cfg_data, edge_file);
    fprintf('Edge table written to: %s\n', edge_file);
end

%% COMMUNICATION MODELS AND DESCRIPTIVE CORRELATIONS
% Communication matrices are computed once and reused by regression and
% spectral workflows. Regression uses transformed/z-scored predictors;
% descriptive correlations and fingerprints use the raw matrices.
if run_communication_models || run_main_regressions || ...
        run_supplementary_single_predictors || run_spectral_analysis
    cm = analysis.compute_communication_models(mats, cfg_data);
else
    cm = [];
end

if run_communication_models
    comm_out = fullfile(output_root, 'communication');
    cfg_comm_tables = cfg_data;
    cfg_comm_tables.out_dir = comm_out;
    analysis.write_communication_tables(cm, mats, cfg_comm_tables);
    save(fullfile(comm_out, 'communication_models.mat'), 'cm', '-v7.3');
end

%% STAGE A NESTED REGRESSION: MAIN COMBINED-MYELIN ANALYSES
% One route-diffusion pair is processed per call. The output writer places
% CSVs under spatial-scale/condition subdirectories. All raw pair results are
% saved for FDR post-processing; correction is outside the two fitting loops.
model_pairs = [
    1 5   % SPE-CMY
    1 6   % SPE-DE
    2 5   % NE-CMY
    2 6   % NE-DE
];

% Columns: display label, interaction, interact_at, effect_mode
conditions = {
    'main-effect',        'none',    'both', 'incremental'
    'interact-caliber',   'caliber', 'both', 'incremental'
    'interact-ed',        'ed',      'both', 'incremental'
    'interact-both',      'both',    'both', 'incremental'
    'total-contribution', 'both',    'both', 'total'
};

if run_main_regressions
    if run_fdr
        assert(cfg_data.use_lower_only, ...
            'Paper FDR requires unique undirected RSN blocks (use_lower_only=true).');
        assert(~isempty(which('bics_fdr_saved_results')) && ...
            ~isempty(which('analysis.write_fdr_tables')), ...
            'Install both FDR helpers and analysis.write_fdr_tables before running.');
        validateattributes(fdr_alpha,{'numeric'},{'scalar','finite','>',0,'<',1});
    end
    condition_keys = {'int_none','int_caliber','int_ed','int_both','total'};
    predictor_tag = strjoin(paper_myelin_predictors,'_');
    main_parent = fullfile(output_root,'regression',['main_results_' predictor_tag]);
    if ~exist(main_parent,'dir'), mkdir(main_parent); end
    [~,run_id] = fileparts(tempname(main_parent));
    main_out = fullfile(main_parent,['run_' datestr(now,'yyyymmdd_HHMMSS') '_' run_id]);
    mkdir(main_out);
    fprintf('\nNEW MAIN-REGRESSION OUTPUT: %s\n',main_out);

    % These four source files collect all conditions in the exact paper order.
    % FDR is applied only after BOTH fitting loops have completed.
    pair_labels = cell(1,size(model_pairs,1));
    paperS = cell(1,size(model_pairs,1));
    for p = 1:size(model_pairs,1)
        pair_labels{p} = strjoin(cfg_data.communication_models(model_pairs(p,:)),'-');
        paperS{p} = struct();
        paperS{p}.meta = struct('label_comm_models',pair_labels{p}, ...
            'CC',model_pairs(p,:),'myelin_predictors',{{paper_myelin_predictors}}, ...
            'base_predictors',{{'binary','caliber'}},'fc_labels',{mats.fc_labels}, ...
            'completed_conditions',{{}},'run_complete',false);
    end

    for c = 1:size(conditions,1)
        for p = 1:size(model_pairs,1)
            cfg_run = cfg_data;
            cfg_run.out_dir = main_out;
            cfg_run.dual_bics_pairs = model_pairs(p,:);
            cfg_run.myelin_predictors = paper_myelin_predictors;
            cfg_run.model_modes = 'all';
            cfg_run.model_levels = [1 2];
            cfg_run.interaction = conditions{c,2};
            cfg_run.interact_at = conditions{c,3};
            cfg_run.effect_mode = conditions{c,4};
            cfg_run.run_stageB = run_stageB;

            pair_label = strjoin(cfg_run.communication_models( ...
                cfg_run.dual_bics_pairs), '-');
            fprintf('\nRegression: %s | %s | combined predictors\n', ...
                conditions{c,1}, pair_label);

            R = analysis.run_nested_regression(cm, mats, pinfo, cfg_run);
            analysis.write_bics_tables(R, cfg_run);

            % Always save raw Stage A outputs, even when Stage B is disabled.
            S = paperS{p};
            S.(condition_keys{c}) = R;
            S.meta.completed_conditions{end+1} = condition_keys{c};
            S.meta.run_complete = c==size(conditions,1);
            S.meta.updated = datestr(now,30);
            paperS{p} = S;
            source_file = fullfile(main_out,['results_total_tmp' pair_labels{p} '.mat']);
            save_paper_source(source_file,S);

            if run_stageB
                mat_dir = fullfile(main_out, 'stageB_mat');
                if ~exist(mat_dir,'dir'), mkdir(mat_dir); end
                save(fullfile(mat_dir, sprintf('%s_%s.mat', ...
                    conditions{c,1}, strrep(pair_label,'-','_'))), ...
                    'R', 'cfg_run', '-v7.3');
            end
        end
    end

    if run_fdr
        F = bics_fdr_saved_results(main_out, ...
            'FilePrefix','results_total_tmp','PairLabels',pair_labels, ...
            'Conditions',condition_keys,'Alpha',fdr_alpha, ...
            'FCLabels',mats.fc_labels,'RSNLabels',pinfo.clabels_short, ...
            'OutputDir',fullfile(main_out,'FDR'),'Save',true,'Overwrite',false);
        analysis.write_fdr_tables(F,fullfile(main_out,'FDR'));
        disp(F.family_summary);
    end
end

%% SUPPLEMENTAL INDIVIDUAL-MYELIN-PREDICTOR ANALYSES
% These are main-effect models with MTsat, g-ratio, and delay tested
% individually. Pairwise and combined predictor modes should be run in
% separate calls rather than mixed within one output file.
if run_supplementary_single_predictors
    supp_out = fullfile(output_root, 'regression', 'supplementary_results');

    for p = 1:size(model_pairs,1)
        cfg_run = cfg_data;
        cfg_run.out_dir = supp_out;
        cfg_run.dual_bics_pairs = model_pairs(p,:);
        cfg_run.myelin_predictors = {'MTsat','gratio','delay'};
        cfg_run.model_modes = 'single';
        cfg_run.model_levels = [1 2];
        cfg_run.interaction = 'none';
        cfg_run.interact_at = 'both';
        cfg_run.effect_mode = 'incremental';
        cfg_run.run_stageB = false;

        pair_label = strjoin(cfg_run.communication_models( ...
            cfg_run.dual_bics_pairs), '-');
        fprintf('\nSupplemental regression: %s | individual predictors\n', ...
            pair_label);

        R = analysis.run_nested_regression(cm, mats, pinfo, cfg_run);
        analysis.write_bics_tables(R, cfg_run);
    end
end

%% COMMUNITY DETECTION: EXACT PAPER PARTITIONS
if run_precomputed_communities
    cfg_cd = cfg.default_community_config(rootdir);
    cfg_cd.out_dir = fullfile(output_root, 'community_detection', 'paper_partitions');
    cfg_cd.networks = {'caliber','MTsat','gratio','delay','rate'};
    cfg_cd.mode = 'precomputed';
    cfg_cd.selection_mode = 'paper';
    cfg_cd.plot_diagnostics = true;
    community_results = analysis.run_community_detection([], cfg_cd); %#ok<NASGU>
end

%% OPTIONAL COMMUNITY-DETECTION RERUN
% This reruns the stochastic coarse-to-fine Louvain workflow. It is not
% required to reproduce the exact bundled paper partitions.
if run_community_rerun
    cfg_cd = cfg.default_community_config(rootdir);
    cfg_cd.out_dir = fullfile(output_root, 'community_detection', 'rerun');
    cfg_cd.networks = community_rerun_networks;
    cfg_cd.mode = 'rerun';
    cfg_cd.selection_mode = 'none';
    cfg_cd.use_parallel = false;
    community_rerun_results = analysis.run_community_detection( ...
        mats, cfg_cd); %#ok<NASGU>
end

%% SPECTRAL ANALYSIS
if run_spectral_analysis
    cfg_spec = cfg.default_spectral_config(rootdir);
    cfg_spec.out_dir = fullfile(output_root, 'spectral_analysis');
    cfg_spec.matched_indices_mode = 'paper';
    spectral_results = analysis.run_spectral_analysis( ...
        mats, cm, pinfo, cfg_spec); %#ok<NASGU>
end

fprintf('\nReplication workflow complete. Outputs: %s\n', output_root);

function save_paper_source(filename,S)
% Replace only the current run's checkpoint after its temporary save succeeds.
tmp = [tempname(fileparts(filename)) '.mat'];
cleanup = onCleanup(@()delete_temp(tmp));
save(tmp,'S','-v7.3');
[ok,msg] = movefile(tmp,filename,'f');
assert(ok,'Replication:SaveFailed','%s',msg);
clear cleanup
end

function delete_temp(path)
if isfile(path), delete(path); end
end
