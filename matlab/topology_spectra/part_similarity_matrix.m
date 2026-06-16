function [AMI, VI] = part_similarity_matrix(C)
% C: N x M matrix of M partitions on the same N nodes (labels as positive integers)
M = size(C,2);
AMI = zeros(M); VI = zeros(M);
for i = 1:M
    for j = i+1:M
        [ami_ij, vi_ij] = ami_vi(C(:,i), C(:,j));
        AMI(i,j) = ami_ij; AMI(j,i) = ami_ij;
        VI(i,j)  = vi_ij;  VI(j,i)  = vi_ij;
    end
end

%--------------------------------------------------------------------------
end