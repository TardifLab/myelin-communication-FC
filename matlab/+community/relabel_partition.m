function partition = relabel_partition(partition, min_size)
%RELABEL_PARTITION Remove small communities and relabel retained groups 1..K.

if nargin < 2 || isempty(min_size), min_size = 1; end
partition = double(partition(:));
labels = unique(partition(isfinite(partition)),'stable');
out = nan(size(partition));
next_label = 0;
for k = 1:numel(labels)
    idx = partition == labels(k);
    if nnz(idx) >= min_size
        next_label = next_label + 1;
        out(idx) = next_label;
    end
end
partition = out;
end
