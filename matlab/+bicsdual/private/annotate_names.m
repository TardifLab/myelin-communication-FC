function [base_name, full_names, red_names, info] = annotate_names(arg1, arg2, varargin)
% ANNOTATE_NAMES v2 — tolerant to 1-arg or 2-arg usage.
%
% Valid forms
%   [base_name, full_names, red_names, info] = annotate_names(base, models, ...)
%   [base_name, full_names, red_names, info] = annotate_names(base)          % base only
%   [base_name, full_names, red_names, info] = annotate_names(models_only)   % no base; models
%
% Name-Value:
%   'ModelLabels', 'Reduced', 'ReducedLabels', 'Prefix'  (as before)

% ---------- Parse name-value ----------
p = inputParser;
p.addParameter('ModelLabels',   {}, @(x) iscell(x) || isstring(x));
p.addParameter('Reduced',       {}, @(x) iscell(x));
p.addParameter('ReducedLabels', {}, @(x) iscell(x) || isstring(x));
p.addParameter('Prefix',        '', @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
mdl_labels   = cellstr(p.Results.ModelLabels);
reduced_list = p.Results.Reduced;
red_labels   = cellstr(p.Results.ReducedLabels);
prefix       = char(p.Results.Prefix);

% ---------- Normalize arguments ----------
base = []; models = [];

if nargin>=2 && ~isempty(arg2)
    % Classic 2-arg form
    base   = arg1;
    models = arg2;
else
    % 1-arg form: try to detect whether it's a base-like struct or models
    if isstruct(arg1) && isfield(arg1,'Xcells')
        base   = arg1;
        models = {};         % base-only
    elseif iscell(arg1)
        % Ambiguous: assume models-only
        base   = struct('Xcells', {{}}, 'labels_base', {{}}); % empty base
        models = arg1;
    else
        error('annotate_names: bad usage. Provide (base, models) or a single base/models argument.');
    end
end

% ---------- Normalize base ----------
if isstruct(base)
    if isfield(base,'Xcells'), BASE = base.Xcells; else, error('base struct missing .Xcells'); end
    if isfield(base,'labels_base') && ~isempty(base.labels_base)
        base_lbls = base.labels_base;
    else
        base_lbls = arrayfun(@(i) sprintf('base%02d',i), 1:numel(BASE), 'uni', 0);
    end
elseif iscell(base)
    BASE     = base;
    base_lbls= arrayfun(@(i) sprintf('base%02d',i), 1:numel(BASE), 'uni', 0);
else
    error('annotate_names: invalid base');
end

% ---------- Sizes ----------
M = numel(models);
if isempty(mdl_labels), mdl_labels = default_model_labels(models); end
hasRED = ~isempty(reduced_list);
if hasRED && numel(reduced_list) ~= M
    error('Reduced list must match number of models.');
end
if ~isempty(red_labels) && numel(red_labels) ~= M
    error('ReducedLabels must match the number of models.');
end

% ---------- Compose names ----------
if isempty(BASE)
    base_name = 'BASE{}';
else
    base_name = sprintf('BASE{%s}', strjoin(base_lbls, ' + '));
end

full_names = cell(1,M);
red_names  = cell(1,M);

info = struct();
info.base_k   = numel(BASE);
info.base_lbl = base_lbls;
info.spec_k   = zeros(1,M);
info.red_k    = zeros(1,M);
info.spec_lbl = cell(1,M);
info.red_lbl  = cell(1,M);

for m = 1:M
    SPECm = models{m};
    info.spec_k(m) = numel(SPECm);

    % SPEC label (prefer provided model label; else synthesize)
    if ~isempty(mdl_labels) && ~isempty(mdl_labels{m})
        spec_str = sanitize_label(mdl_labels{m});
    else
        spec_str = synthesize_label('SPEC', SPECm);
    end
    info.spec_lbl{m} = spec_str;

    if ~isempty(prefix)
        full_names{m} = sprintf('%s: SPEC{%s}', prefix, spec_str);
    else
        full_names{m} = sprintf('SPEC{%s}', spec_str);
    end

    % RED label if present
    if hasRED && ~isempty(reduced_list{m})
        REDm = reduced_list{m};
        info.red_k(m) = numel(REDm);

        if ~isempty(red_labels) && ~isempty(red_labels{m})
            red_str = sanitize_label(red_labels{m});
        else
            red_str = synthesize_label('RED', REDm);
        end
        info.red_lbl{m} = red_str;

        if ~isempty(prefix)
            red_names{m} = sprintf('%s: RED{%s}', prefix, red_str);
        else
            red_names{m} = sprintf('RED{%s}', red_str);
        end
    else
        red_names{m} = '';
    end
end

% ---------- local helpers ----------
function L = default_model_labels(MDL)
    L = cell(1,numel(MDL));
    for ii = 1:numel(MDL)
        K = numel(MDL{ii});
        L{ii} = sprintf('K=%d', K);
    end
end

function s = synthesize_label(tag, CELL)
    K = numel(CELL);
    s = sprintf('%s(K=%d)', tag, K);
end

function s = sanitize_label(s0)
    s = strtrim(char(s0));
    s = strrep(s, '{', '');
    s = strrep(s, '}', '');
    s = regexprep(s, '\s+', ' ');
end
end
