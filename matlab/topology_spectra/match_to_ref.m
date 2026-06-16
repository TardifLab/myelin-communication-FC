function out = match_to_ref(Yref, Cpart)
% MATCH_TO_REF  Contingency vs reference + Hungarian match + safe column order.
% Works with MUNKRES versions that return either a vector (1xR) or a matrix (R x K).
%
% Inputs
%   Yref   : N x 1 reference labels (positive ints; NaN/0 ignored)
%   Cpart  : N x 1 partition labels  (positive ints; NaN/0 ignored)
%
% Outputs (struct)
%   .CM           : R x K contingency (counts)
%   .Cmap_row     : R x K row-normalized contingency (rows sum to 1)
%   .Cmap_col     : R x K column-normalized contingency (cols sum to 1)
%   .row_labels   : 1 x R original Yref IDs (e.g., [1 2 3 4 5 6 7])
%   .col_labels   : 1 x K original Cpart community IDs
%   .row_sizes    : R x 1 sizes of reference communities
%   .col_sizes    : 1 x K sizes of target communities
%   .assign_rc    : 1 x R row->col vector (0 = unassigned)
%   .col_order    : 1 x K permutation of columns for plotting (assigned first)
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------


% -------- sanitize labels (ignore NaN/0) --------
Yref  = Yref(:);  Cpart = Cpart(:);
valid = ~isnan(Yref) & Yref>0 & ~isnan(Cpart) & Cpart>0;
Yref  = Yref(valid);
Cpart = Cpart(valid);

% -------- compact relabel, preserve originals --------
[Yr_unique, ~, Yr_idx] = unique(Yref,'stable');   % R rows (original IDs)
[Cp_unique, ~, Cp_idx] = unique(Cpart,'stable');  % K cols (original IDs)
R = numel(Yr_unique);
K = numel(Cp_unique);

% -------- contingency: rows=Yref groups, cols=partition groups --------
CM = accumarray([Yr_idx, Cp_idx], 1, [R, K]);

% sizes
row_sizes = sum(CM, 2);           % R x 1
col_sizes = sum(CM, 1);           % 1 x K

% -------- normalized maps --------
rs = row_sizes; rs(rs==0) = 1;    % guard divide-by-zero
cs = col_sizes; cs(cs==0) = 1;

Cmap_row = CM ./ rs;              % each row sums to 1 (coverage of Yref by Cpart)
Cmap_col = CM ./ cs;              % each col sums to 1 (purity of Cpart wrt Yref)

% -------- Hungarian assignment (maximize overlap) --------
cost = max(CM(:)) - CM;           % convert to min-cost
assign_any = munkres(cost);       % your munkres returns matrix (R x K) logical

% normalize to 1xR row->col vector
if isvector(assign_any)
    assign_rc = assign_any(:).';                % 1 x R
    if numel(assign_rc) ~= R
        error('munkres returned vector length %d; expected %d (R).', numel(assign_rc), R);
    end
else
    if ~isequal(size(assign_any), [R K])
        error('munkres returned %s; expected [%d %d] or [1 %d].', ...
              mat2str(size(assign_any)), R, K, R);
    end
    A = logical(assign_any);
    assign_rc = zeros(1,R);
    for r = 1:R
        c = find(A(r,:), 1, 'first');          % 0 if row unassigned
        if isempty(c), c = 0; end
        assign_rc(r) = c;
    end
end

% -------- build a proper 1xK column permutation for plotting --------
assigned_cols = assign_rc(assign_rc > 0);
assigned_cols = unique(assigned_cols,'stable');       % keep first occurrences
remaining_cols = setdiff(1:K, assigned_cols, 'stable');
col_order = [assigned_cols(:).' remaining_cols(:).']; % all columns exactly once

% -------- package --------
out = struct();
out.CM         = CM;
out.Cmap_row   = Cmap_row;
out.Cmap_col   = Cmap_col;
out.row_labels = Yr_unique(:).';
out.col_labels = Cp_unique(:).';
out.row_sizes  = row_sizes;
out.col_sizes  = col_sizes;
out.assign_rc  = assign_rc;
out.col_order  = col_order;

%--------------------------------------------------------------------------
end
