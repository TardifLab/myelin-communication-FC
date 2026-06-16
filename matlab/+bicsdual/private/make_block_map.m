function [blk_ij, blk_labels] = make_block_map(Nntwk, use_lower, pinfo)
ij = [];
for i = 1:Nntwk
    if use_lower, jset = 1:i; else, jset = 1:Nntwk; end
    for j = jset, ij = [ij; i j]; end %#ok<AGROW>
end
blk_ij = ij;
blk_labels = cell(size(blk_ij,1),1);
for k = 1:size(blk_ij,1)
    i = blk_ij(k,1); j = blk_ij(k,2);
    blk_labels{k} = sprintf('%s–%s', label_of(i,pinfo), label_of(j,pinfo));
end
end
