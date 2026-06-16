function sa_axis = resolve_sa_axis(config, pinfo, N)
%RESOLVE_SA_AXIS Resolve the sensory-association axis for spectral analyses.
%
% Priority:
%   1. config.sa_axis
%   2. compatible field in pinfo
%   3. numeric vector loaded from config.sa_axis_file

sa_axis = [];

%% 1. Explicit vector supplied in config
if isfield(config, 'sa_axis') && ~isempty(config.sa_axis)
    sa_axis = config.sa_axis;
end

%% 2. Try node metadata
if isempty(sa_axis) && ~isempty(pinfo)

    candidate_fields = { ...
        'sa_axis', ...
        'SA_axis', ...
        'sensory_association', ...
        'sensory_association_axis', ...
        'SA'};

    if istable(pinfo)
        vars = pinfo.Properties.VariableNames;

        for k = 1:numel(candidate_fields)
            if ismember(candidate_fields{k}, vars)
                sa_axis = pinfo.(candidate_fields{k});
                break
            end
        end

    elseif isstruct(pinfo)
        for k = 1:numel(candidate_fields)
            if isfield(pinfo, candidate_fields{k})
                sa_axis = pinfo.(candidate_fields{k});
                break
            end
        end
    end
end

%% 3. Load vector from .mat file
if isempty(sa_axis) && ...
        isfield(config, 'sa_axis_file') && ...
        ~isempty(config.sa_axis_file)

    if ~exist(config.sa_axis_file, 'file')
        error('SA-axis file not found: %s', config.sa_axis_file);
    end

    S = load(config.sa_axis_file);
    names = fieldnames(S);

    % Prefer clearly named variables
    preferred_names = { ...
        'sa_axis', ...
        'SA_axis', ...
        'sensory_association', ...
        'sensory_association_axis', ...
        'SA'};

    for k = 1:numel(preferred_names)
        if isfield(S, preferred_names{k})
            candidate = S.(preferred_names{k});

            if isnumeric(candidate) && isvector(candidate)
                sa_axis = candidate;
                break
            end
        end
    end

    % Otherwise use the first numeric vector of the expected length
    if isempty(sa_axis)
        for k = 1:numel(names)
            candidate = S.(names{k});

            if isnumeric(candidate) && ...
                    isvector(candidate) && ...
                    numel(candidate) == N
                sa_axis = candidate;
                break
            end
        end
    end

    if isempty(sa_axis)
        error(['No numeric vector of length %d was found in:\n%s\n', ...
               'Variables present: %s'], ...
               N, config.sa_axis_file, strjoin(names, ', '));
    end
end

%% Validate
if isempty(sa_axis)
    error(['No sensory-association axis was supplied. Set ', ...
           'config.sa_axis or config.sa_axis_file.']);
end

sa_axis = double(sa_axis(:));

if numel(sa_axis) ~= N
    error('SA axis has %d values; expected %d.', numel(sa_axis), N);
end

if any(~isfinite(sa_axis))
    error('SA axis contains NaN or Inf values.');
end

end
