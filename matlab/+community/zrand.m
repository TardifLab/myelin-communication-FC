function z = zrand(a, b)
%ZRAND Compute the z-score of the Rand index between two partitions.
%
% Community labels are converted to double because some versions of
% community_louvain return integer arrays that are incompatible with
% dummyvar inside fcn_randz.

a = double(a(:));
b = double(b(:));

% Remove rows that are invalid in either partition.
ok = isfinite(a) & isfinite(b);
a = a(ok);
b = b(ok);

if isempty(a) || numel(a) ~= numel(b)
    z = NaN;
    return
end

% Relabel communities consecutively from 1...K.
% This avoids problems with zero, skipped, or nonconsecutive labels.
[~, ~, a] = unique(a, 'stable');
[~, ~, b] = unique(b, 'stable');

a = double(a);
b = double(b);

z = fcn_randz(a, b);

end
