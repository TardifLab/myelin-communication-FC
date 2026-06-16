function [Xlist, labels] = predictor_pack(U, idx, predictors, label_prefix)
% PREDICTOR_PACK  General helper to gather arbitrary predictors for a
%                 given route/diff index (generalizing myelin_pack).
%
%   [Xlist, labels] = predictor_pack(U, idx, predictors)
%   [Xlist, labels] = predictor_pack(U, idx, predictors, label_prefix)
%
% Inputs
%   U          : struct with fields named by predictor type
%                (e.g., 'MTsat', 'gratio', 'delay', 'caliber', ...).
%                Each field is expected to be a 1xN cell array where
%                U.(field){idx} is an [Nedge x K] matrix of predictors.
%
%   idx        : scalar index selecting which route/diff entry to use
%                (e.g., iR or iD in the dual BICS framework).
%
%   predictors : cellstr or char specifying which fields in U to pack.
%                Examples:
%                   {'MTsat','gratio'}
%                   {'caliber'}
%                This is intended to replace / generalize myelin_pack,
%                allowing use for myelin metrics & caliber
%
%   label_prefix (optional) : char or string to prepend to each label.
%                If provided and non-empty, labels will be:
%                   '<label_prefix>_<predictor>'
%                Otherwise labels are just the predictor names.
%
% Outputs
%   Xlist      : 1 x P cell of predictor matrices, where P is
%                numel(predictors). Each Xlist{p} is taken from
%                U.(predictors{p}){idx}.
%
%   labels     : 1 x P cell of char labels corresponding to predictors
%                (with optional prefix).
%
% Notes
%   - This is a generalization of your myelin_pack: instead of being
%     restricted to myelin fields, it will happily pack caliber as well,
%     or any other field in U that follows the U.field{idx} convention.
%   - You can use this both for "myelin stacks" and for caliber-only
%     stacks by simply changing the predictors list.

    % Optional prefix
    if nargin < 4
        label_prefix = '';
    end
    if isstring(label_prefix)
        label_prefix = char(label_prefix);
    end

    % Normalize predictor list to a cellstr
    if ischar(predictors) || isstring(predictors)
        predictors = cellstr(predictors);
    elseif ~iscell(predictors)
        error('predictor_pack:BadPredictors',...
              'predictors must be a char, string, or cellstr.');
    end

    nP = numel(predictors);
    Xlist  = cell(1, nP);
    labels = cell(1, nP);

    % Basic index sanity
    if ~isscalar(idx) || idx ~= floor(idx) || idx < 1
        error('predictor_pack:BadIndex',...
              'idx must be a positive scalar integer.');
    end

    for p = 1:nP
        name = predictors{p};

        % Check field exists
        if ~isfield(U, name)
            error('predictor_pack:MissingField',...
                  'U is missing predictor field ''%s''.', name);
        end

        fld = U.(name);

        % Expect a cell array (U.field{idx})
        if ~iscell(fld)
            error('predictor_pack:NotCell',...
                  'U.%s must be a cell array; got %s.', ...
                  name, class(fld));
        end

        if idx > numel(fld)
            error('predictor_pack:IndexOutOfRange',...
                  'idx=%d is out of range for U.%s (numel=%d).', ...
                  idx, name, numel(fld));
        end

        Xi = fld{idx};

        if isempty(Xi)
            warning('predictor_pack:EmptyPredictor',...
                    'U.%s{%d} is empty.', name, idx);
        end

        Xlist{p} = Xi;

        if isempty(label_prefix)
            labels{p} = name;
        else
            labels{p} = sprintf('%s_%s', label_prefix, name);
        end
    end
end
