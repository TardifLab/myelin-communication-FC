function varargout = sanitize_naninf_rows(varargin)
% SANITIZE_NANINF_ROWS  Drop rows containing NaN/Inf across multiple arrays.
%
% Usage:
%   [Xb_c, Xf_c, y_c, info] = sanitize_naninf_rows(Xb, Xf, y);
%   [...] = sanitize_naninf_rows(X1, X2, ..., 'dim', 1);  % observations in rows (default)
%   [...] = sanitize_naninf_rows(X1, X2, ..., 'dim', 2);  % observations in columns
%
% Inputs:
%   - Any number of numeric/sparse arrays with the same # of observations along DIM
%   - Name-value:
%       'dim' : 1 (default) or 2  — along which dimension observations lie
%
% Outputs:
%   - Cleaned arrays in the same order (rows/cols with any NaN/Inf dropped)
%   - Final output is an info struct with fields:
%       .dim, .n_obs_in, .n_obs_out, .n_removed, .mask_keep, .mask_drop
%       .per_input.nan_count, .per_input.inf_count
%
% Notes:
%   - This helper is safe for dense/sparse inputs.
%   - If all rows are dropped, returns empty arrays with 0 observations.

    % -------- parse name-value (only 'dim') --------
    dim = 1;
    nv_start = numel(varargin);
    while nv_start >= 2 && (ischar(varargin{nv_start-1}) || isStringScalar(varargin{nv_start-1}))
        key = lower(string(varargin{nv_start-1}));
        val = varargin{nv_start};
        if key == "dim"
            if ~(isscalar(val) && (val==1 || val==2))
                error('sanitize_naninf_rows:BadDim','"dim" must be 1 or 2.');
            end
            dim = val;
            nv_start = nv_start - 2;
        else
            error('sanitize_naninf_rows:UnknownNV','Unknown parameter "%s".', key);
        end
    end
    n_inputs = nv_start;

    if n_inputs < 1
        error('sanitize_naninf_rows:NoInputs','Provide at least one array.');
    end

    % -------- basic checks & infer n_obs --------
    sz = cell(n_inputs,1);
    for i = 1:n_inputs
        if ~(isnumeric(varargin{i}) || issparse(varargin{i}))
            error('sanitize_naninf_rows:Type','All inputs must be numeric/sparse.');
        end
        sz{i} = size(varargin{i});
    end

    % Determine # observations along dim; ensure consistent
    n_obs = sz{1}(dim);
    for i = 2:n_inputs
        if sz{i}(dim) ~= n_obs
            error('sanitize_naninf_rows:SizeMismatch', ...
                  'Input %d has %d obs along dim=%d; expected %d.', ...
                  i, sz{i}(dim), dim, n_obs);
        end
    end

    % -------- build keep mask: finite across ALL inputs --------
    mask_keep = true(n_obs,1);
    per_input.nan_count = zeros(n_inputs,1);
    per_input.inf_count = zeros(n_inputs,1);

    for i = 1:n_inputs
        A = varargin{i};
        % Logical mask of finite values per observation
        if dim == 1
            % rows are observations
            finite_row = all(isfinite(full(A)), 2);
        else
            % columns are observations
            finite_row = all(isfinite(full(A)), 1).';
        end
        % Update per-input counts
        if dim == 1
            per_input.nan_count(i) = sum(any(isnan(full(A)), 2));
            per_input.inf_count(i) = sum(any(isinf(full(A)), 2));
        else
            per_input.nan_count(i) = sum(any(isnan(full(A)), 1));
            per_input.inf_count(i) = sum(any(isinf(full(A)), 1));
        end
        % Accumulate across inputs (keep only obs that are finite in ALL)
        mask_keep = mask_keep & finite_row;
    end

    % -------- apply mask to each input --------
    varargout = cell(1, n_inputs + 1);
    for i = 1:n_inputs
        A = varargin{i};
        if dim == 1
            varargout{i} = A(mask_keep, :);
        else
            varargout{i} = A(:, mask_keep);
        end
    end

    % -------- info struct --------
    info = struct();
    info.dim        = dim;
    info.n_obs_in   = n_obs;
    info.n_obs_out  = sum(mask_keep);
    info.n_removed  = info.n_obs_in - info.n_obs_out;
    info.mask_keep  = mask_keep;
    info.mask_drop  = ~mask_keep;
    info.per_input  = per_input;

    varargout{n_inputs+1} = info;
end
