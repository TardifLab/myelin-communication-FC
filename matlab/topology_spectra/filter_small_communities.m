function [C_out,info] = filter_small_communities(C_in,min_size)
% FILTER_SMALL_COMMUNITIES  Remove tiny communities and renumber.
% 
% Inputs
%   C_in      : N x P matrix of integer community labels (NaN/0 allowed)
%   min_size  : scalar (e.g., 2). Communities with < min_size nodes are removed (set to NaN).
%
% Outputs
%   C_out     : N x P filtered/renumbered labels (tiny comms -> NaN; remaining -> 1..k)
%   info      : struct with per-partition details:
%                 .removed_labels   : labels removed due to size
%                 .removed_sizes    : their sizes
%                 .kept_labels_old  : surviving labels (original IDs)
%                 .kept_labels_new  : their new IDs (1..k after filtering)
%                 .k_new            : final # communities after filtering
%
% Notes
% - NaN entries are left as NaN.
% - Label 0 (if present) is treated as background and preserved as 0 (not renumbered).
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------


if nargin < 2 || isempty(min_size), min_size = 2; end

[N, P] = size(C_in);
C_out  = C_in;

% preallocate info
info = repmat(struct('removed_labels',[], 'removed_sizes',[], ...
                     'kept_labels_old',[], 'kept_labels_new',[], ...
                     'k_new',0), P, 1);

for j = 1:P
    x = C_out(:,j);
    valid = ~isnan(x) & x > 0;                % positive labels only (ignore 0/NaN)

    if ~any(valid)
        % nothing to filter or renumber
        info(j).k_new = 0;
        continue
    end

    % counts per positive label
    labs = unique(x(valid), 'stable');
    cnts = arrayfun(@(L) sum(x==L), labs);

    % tiny communities -> NaN
    tiny_mask = cnts < min_size;
    tiny_labs = labs(tiny_mask);
    if ~isempty(tiny_labs)
        x(ismember(x, tiny_labs)) = NaN;
    end

    % survivors: renumber to 1..k in stable label order
    valid2 = ~isnan(x) & x > 0;
    if any(valid2)
        labs2 = unique(x(valid2),'stable');      % surviving labels (old IDs)
        k_new = numel(labs2);

        % map old -> new using ismember indexing (fast, no containers.Map needed)
        [tf, loc] = ismember(x, labs2);          % loc gives position in labs2 or 0
        % preserve zeros and NaNs
        x_new = NaN(size(x));
        x_new(x==0) = 0;                          % keep background zeros as 0
        idx = tf & x>0;                           % survivors (positive labels)
        x_new(idx) = loc(idx);                    % assign 1..k

        x = x_new;                                % replace
    else
        labs2 = [];
        k_new = 0;
        % keep 0s and NaNs as they are; all positive labels were removed
        x(~isnan(x) & x>0) = NaN;
    end

    % write back
    C_out(:,j) = x;

    % bookkeeping
    info(j).removed_labels  = tiny_labs(:);
    info(j).removed_sizes   = cnts(tiny_mask);
    info(j).kept_labels_old = labs2(:);
    if isempty(labs2)
        info(j).kept_labels_new = [];
    else
        info(j).kept_labels_new = (1:numel(labs2)).';
    end
    info(j).k_new = k_new;
end

%--------------------------------------------------------------------------
end
