function [Xlist, labels] = base_pack(U, idx, base_predictors, label_suffix)
% BASE_PACK  Gather base predictors for a given route/diff index.
%
%   [Xlist, labels] = base_pack(U, idx, base_predictors)
%
% Inputs
%   U              : struct with fields named by predictor type
%                    (e.g., 'binary', 'caliber', 'MTsat', 'gratio', ...).
%                    Each field is expected to be a 1xN cell array where
%                    U.(field){idx} is an [Nedge x K] matrix of predictors.
%
%   idx            : scalar index selecting which route/diff entry to use
%                    (e.g., iR or iD in the dual BICS framework).
%
%   base_predictors: cellstr or char specifying which fields in U to
%                    include in the base. Example:
%                       {'binary','caliber'}
%
%   label_suffix (optional) : char or string to append to each label.
%                If provided and non-empty, labels will be:
%                   '<predictor>_<label_prefix>'
%                Otherwise labels are just the predictor names.
%
% Outputs
%   Xlist          : 1 x P cell of predictor matrices, where P is
%                    numel(base_predictors). Each Xlist{p} is taken
%                    from U.(base_predictors{p}){idx}.
%
%   labels         : 1 x P cell of char labels mirroring base_predictors.
%
% Notes
%   - ED is *not* handled here. EDis appended separately
%     in the orchestration code.

    % Optional suffix
    if nargin < 4
        label_suffix = '';
    end
    if isstring(label_suffix)
        label_suffix = char(label_suffix);
    end

    % Normalize predictor list to a cellstr
    if ischar(base_predictors) || isstring(base_predictors)
        base_predictors = cellstr(base_predictors);
    elseif ~iscell(base_predictors)
        error('base_pack:BadPredictors',...
              'base_predictors must be a char, string, or cellstr.');
    end

    nP = numel(base_predictors);
    Xlist  = cell(1, nP);
    labels = cell(1, nP);

    % Basic index sanity
    if ~isscalar(idx) || idx ~= floor(idx) || idx < 1
        error('base_pack:BadIndex',...
              'idx must be a positive scalar integer.');
    end

    for p = 1:nP
        name = base_predictors{p};

        % Check field exists
        if ~isfield(U, name)
            error('base_pack:MissingField',...
                  'U is missing predictor field ''%s''.', name);
        end

        fld = U.(name);

        % Expect a cell array (U.field{idx})
        if ~iscell(fld)
            error('base_pack:NotCell',...
                  'U.%s must be a cell array; got %s.', ...
                  name, class(fld));
        end

        if idx > numel(fld)
            error('base_pack:IndexOutOfRange',...
                  'idx=%d is out of range for U.%s (numel=%d).', ...
                  idx, name, numel(fld));
        end

        Xi = fld{idx};

        if isempty(Xi)
            warning('base_pack:EmptyPredictor',...
                    'U.%s{%d} is empty.', name, idx);
        end

        Xlist{p}  = Xi;

        if isempty(label_suffix)
            labels{p} = name;
        else
            labels{p} = sprintf('%s_%s', name, label_suffix);
        end
    end
end
