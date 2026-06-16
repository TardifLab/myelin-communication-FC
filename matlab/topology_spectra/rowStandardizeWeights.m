function W = rowStandardizeWeights(D)
% D: NxN weighted connectivity (can be asymmetric). Missing = NaN or 0.
% - each non-empty row of W sums to 1, diagonal is 0.
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------


    % 1) Replace NaNs with 0 and zero the diagonal
    D = D;
    D(isnan(D)) = 0;
    D(1:size(D,1)+1:end) = 0;

    % (Optional) Symmetrize
    % D = max(D, D.');   % or D = (D + D.')/2;

    % 2) Row sums
    rs = sum(D, 2);

    % 3) Row-standardize (leave isolated rows as all zeros)
    W = zeros(size(D), 'like', D);
    nz = rs > 0;
    W(nz, :) = D(nz, :) ./ rs(nz);

end
