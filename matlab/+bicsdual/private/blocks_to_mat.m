function M = blocks_to_mat(bvec, blk_ij, Nntwk, symmetrize)
M = zeros(Nntwk,Nntwk);
for b=1:size(blk_ij,1)
    i=blk_ij(b,1); j=blk_ij(b,2);
    M(i,j) = bvec(b);
    if symmetrize, M(j,i) = bvec(b); end
end
end