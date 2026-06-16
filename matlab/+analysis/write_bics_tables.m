function write_bics_tables(results, cfg)
%WRITE_BICS_TABLES Write BICS results to organized CSV files.
%
% Output hierarchy:
%   cfg.out_dir/
%       global/<condition>/bics_deltaR2_<pair>_<mode>_<levels>.csv
%       network/<condition>/bics_deltaR2_<pair>_<mode>_<levels>.csv
%       node/<condition>/bics_deltaR2_<pair>_<mode>_<levels>.csv
%
% where <condition> is one of:
%   main-effect
%   interact-caliber
%   interact-ed
%   interact-both
%   total-contribution

if ~exist(cfg.out_dir,'dir')
    mkdir(cfg.out_dir);
end

T = analysis.results_to_tables(results.stats_all, results.fc_labels);

condition_label = local_condition_label(cfg);
file_suffix = local_make_file_suffix(results, cfg);

if isfield(T,'global') && ~isempty(T.global)
    T.global = local_add_metadata(T.global, results, cfg, condition_label);
    local_write_table(T.global, cfg.out_dir, 'global', condition_label, file_suffix);
end

if isfield(T,'network') && ~isempty(T.network)
    T.network = local_add_metadata(T.network, results, cfg, condition_label);
    local_write_table(T.network, cfg.out_dir, 'network', condition_label, file_suffix);
end

if isfield(T,'node') && ~isempty(T.node)
    T.node = local_add_metadata(T.node, results, cfg, condition_label);
    local_write_table(T.node, cfg.out_dir, 'node', condition_label, file_suffix);
end

end


function local_write_table(T, out_root, scale_label, condition_label, file_suffix)

outdir = fullfile(out_root, scale_label, condition_label);

if ~exist(outdir, 'dir')
    mkdir(outdir);
end

outfile = fullfile(outdir, ['bics_deltaR2' file_suffix '.csv']);
writetable(T, outfile);

end


function suffix = local_make_file_suffix(results, cfg)

labels = {};

if isfield(results, 'communication_models') && ~isempty(results.communication_models)
    labels = results.communication_models;
elseif isfield(cfg, 'dual_bics_pairs') && isfield(cfg, 'communication_models')
    pair = cfg.dual_bics_pairs(1,:);
    labels = cfg.communication_models(pair);
end

if ischar(labels) || isstring(labels)
    labels = cellstr(labels);
end

suffix_parts = {};

% Communication pair, e.g. SPE_CMY
if ~isempty(labels)
    suffix_parts{end+1} = strjoin(labels, '_'); %#ok<AGROW>
end

% Myelin predictor grouping, e.g. all / single / pairs
if isfield(cfg, 'model_modes')
    mode_label = cfg.model_modes;
    if iscell(mode_label)
        mode_label = strjoin(mode_label, '_');
    end
    suffix_parts{end+1} = char(string(mode_label)); %#ok<AGROW>
end

% Model levels, e.g. levels-L1L2
if isfield(cfg, 'model_levels')
    lvl = cfg.model_levels;

    if isnumeric(lvl)
        if isequal(lvl, [1 2])
            lvl_label = 'levels-L1L2';
        elseif isequal(lvl, 1)
            lvl_label = 'levels-L1';
        elseif isequal(lvl, 2)
            lvl_label = 'levels-L2';
        else
            lvl_label = ['levels-L' strjoin(cellstr(string(lvl)), '_')];
        end
    else
        lvl_label = ['levels-' char(string(lvl))];
    end

    suffix_parts{end+1} = lvl_label; %#ok<AGROW>
end

if isempty(suffix_parts)
    suffix = '';
else
    suffix = ['_' strjoin(suffix_parts, '_')];
end

% Sanitize for filenames
suffix = regexprep(suffix, '[^A-Za-z0-9_]+', '_');

end


function T = local_add_metadata(T, results, cfg, condition_label)

n = height(T);

route_model = strings(n,1);
diffusion_model = strings(n,1);
communication_pair = strings(n,1);
condition = repmat(string(condition_label), n, 1);

labels = {};

if isfield(results, 'communication_models') && ~isempty(results.communication_models)
    labels = results.communication_models;
elseif isfield(cfg, 'dual_bics_pairs') && isfield(cfg, 'communication_models')
    pair = cfg.dual_bics_pairs(1,:);
    labels = cfg.communication_models(pair);
end

if ischar(labels) || isstring(labels)
    labels = cellstr(labels);
end

if numel(labels) >= 1
    route_model(:) = string(labels{1});
end

if numel(labels) >= 2
    diffusion_model(:) = string(labels{2});
end

if ~isempty(labels)
    communication_pair(:) = string(strjoin(labels, '_'));
end

T = addvars(T, route_model, diffusion_model, communication_pair, condition, ...
    'Before', 1);

end


function condition_label = local_condition_label(cfg)

% Default to incremental main-effect model if fields are absent
effect_mode = 'incremental';
interaction = 'none';

if isfield(cfg, 'effect_mode') && ~isempty(cfg.effect_mode)
    effect_mode = lower(char(string(cfg.effect_mode)));
end

if isfield(cfg, 'interaction') && ~isempty(cfg.interaction)
    interaction = lower(char(string(cfg.interaction)));
end

% Total contribution takes precedence because it describes the effect being
% tested, regardless of the interaction setting.
if strcmp(effect_mode, 'total')
    condition_label = 'total-contribution';
    return
end

switch interaction
    case {'none', 'main', 'main_effect'}
        condition_label = 'main-effect';

    case {'caliber', 'calibre'}
        condition_label = 'interact-caliber';

    case {'ed', 'euclidean', 'distance', 'euclidean_distance'}
        condition_label = 'interact-ed';

    case {'both', 'caliber_ed', 'ed_caliber'}
        condition_label = 'interact-both';

    otherwise
        condition_label = ['interact-' interaction];
end

end
