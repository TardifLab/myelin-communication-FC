function results = run_nested_regression(cm, mats, pinfo, cfg)
%RUN_NESTED_REGRESSION Run nested regression.

% Check model modes
if iscell(cfg.model_modes) && numel(cfg.model_modes) > 1
    error(['Please use one cfg.model_modes value per run: ', ...
        '''single'', ''pairs'', or ''all''. ', ...
        'Run separate configs for main and supplemental analyses.']);
end

% Check interaction/effect settings
valid_interactions = {'none','caliber','ed','both'};
valid_effect_modes = {'incremental','total'};

if ~ismember(lower(cfg.interaction), valid_interactions)
    error(['Invalid cfg.interaction = "%s". Use one of: %s. ', ...
           'For total contribution, set cfg.effect_mode = ''total'', ', ...
           'not cfg.interaction = ''total''.'], ...
           cfg.interaction, strjoin(valid_interactions, ', '));
end

if ~ismember(lower(cfg.effect_mode), valid_effect_modes)
    error('Invalid cfg.effect_mode = "%s". Use one of: %s.', ...
          cfg.effect_mode, strjoin(valid_effect_modes, ', '));
end

% Begin
ALLY = cell(size(mats.FC));
for f = 1:numel(mats.FC)
    ALLY{f} = nzzscore(mats.FC{f}, 0);
end

if ~isfield(pinfo,'clabels_short')
    pinfo.clabels_short = pinfo.clabels;
end

tests = struct('myelin_predictors', {cfg.myelin_predictors}, ...
               'label_comm_models', strjoin(cfg.communication_models(cfg.dual_bics_pairs), '-'), ...
               'modes', cfg.model_modes, ...
               'levels', cfg.model_levels, ...
               'interaction', cfg.interaction, ...
               'interact_at', cfg.interact_at, ...
               'effect_mode', cfg.effect_mode);
optsA = struct('use_lower_only', cfg.use_lower_only, 'alpha', cfg.alpha);

[Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, cm.ALLX, cfg.dual_bics_pairs, tests, pinfo, optsA);

results = struct();
results.Delta_all = Delta_all;
results.P_all = P_all;
results.stats_all = stats_all;
results.tests = tests;
results.fc_labels = mats.fc_labels;

if isfield(cfg,'run_stageB') && cfg.run_stageB
    optsB = cfg.stageB;
    optsB.use_lower_only = cfg.use_lower_only;
    [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, cm.ALLX, cfg.dual_bics_pairs, tests, pinfo, Delta_all, optsB);
    results.SI_all = SI_all;
    results.elasticity_all = elasticity_all;
    results.pert_all = pert_all;
end
end
