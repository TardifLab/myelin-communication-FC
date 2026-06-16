function inputs = get_inputs(mats, cm, config)
%GET_INPUTS Assemble structural and communication inputs in paper order.

required_mats = {'caliber','MTsat','gratio'};
for i = 1:numel(required_mats)
    if ~isfield(mats, required_mats{i})
        error('mats.%s is required for spectral analysis.', required_mats{i});
    end
end
if ~isfield(cm, 'raw')
    error('cm.raw is required. Run analysis.compute_communication_models first.');
end

if isfield(cm, 'normalized_inputs')
    normalized = cm.normalized_inputs;
else
    normalized = struct();
end

inputs = struct();
inputs.structural = struct();
inputs.structural.caliber = local_get_structural(normalized, mats, 'caliber');
inputs.structural.MTsat   = local_get_structural(normalized, mats, 'MTsat');
inputs.structural.gratio  = local_get_structural(normalized, mats, 'gratio');

% The manuscript's "delay" structural comparison used the signaling-rate
% matrix derived from tract length and g-ratio.
if isfield(normalized, 'delay_rate')
    inputs.structural.delay = local_clean(normalized.delay_rate);
elseif isfield(mats, 'delay_rate')
    inputs.structural.delay = local_clean(mats.delay_rate);
elseif isfield(mats, 'rate')
    inputs.structural.delay = local_clean(mats.rate);
else
    if ~isfield(mats,'length') || ~isfield(mats,'gratio')
        error(['Delay-derived rate network unavailable: provide ', ...
            'cm.normalized_inputs.delay_rate, mats.delay_rate, mats.rate, ', ...
            'or mats.length plus mats.gratio.']);
    end
    if exist('prep.make_delay','file') ~= 2
        error('prep.make_delay is required to reconstruct the delay/rate network.');
    end
    [~, rate] = prep.make_delay(mats.length, mats.gratio, config.delay);
    inputs.structural.delay = local_clean(rate);
end

% spectral_fingerprint_windows expects 1 x J cells, each containing 1 x M
% communication matrices. Raw matrices reproduce the manuscript workflow.
names = config.dataset_labels;
inputs.communication = cell(1, numel(names));
for j = 1:numel(names)
    name = names{j};
    if ~isfield(cm.raw, name)
        error('cm.raw.%s is missing.', name);
    end
    inputs.communication{j} = cm.raw.(name);
end

if isfield(cm,'model_labels') && ~isempty(cm.model_labels)
    inputs.model_labels = cm.model_labels;
else
    inputs.model_labels = config.communication_models;
end

if numel(inputs.model_labels) ~= numel(config.operator_for_model)
    error('operator_for_model must contain one entry per communication model.');
end
if numel(inputs.communication) ~= numel(config.dataset_labels)
    error('dataset_labels must contain one label per structural dataset.');
end

N = size(inputs.structural.caliber,1);
all_names = fieldnames(inputs.structural);
for i = 1:numel(all_names)
    X = inputs.structural.(all_names{i});
    if ~isequal(size(X), [N N])
        error('Structural matrix %s is not %d x %d.', all_names{i}, N, N);
    end
end
end

function X = local_get_structural(normalized, mats, name)
if isfield(normalized, name)
    X = normalized.(name);
else
    X = mats.(name);
    if exist('prep.normalize_weights','file') == 2
        X = prep.normalize_weights(X);
    end
end
X = local_clean(X);
end

function X = local_clean(X)
X = double(X);
X(~isfinite(X)) = 0;
X(X < 0) = 0;
X = (X + X.') ./ 2;
X(1:size(X,1)+1:end) = 0;
end
