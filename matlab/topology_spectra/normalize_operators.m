function [Ahat,Lsym] = normalize_operators(A)
% Outputs normalized versions of A & L for spectral analysis
% 
% Input:
%   A       : weighted adjacency matrix (NxN)
%
% Output: 
%   Ahat    : symmetric degree-normalized adjacency (spectrum [-1 1])
%   Lsym    : normalized Laplacian (eigenvalues in range [0 2]
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

    % Symmetrize
    A = (A + A')/2;

    % Degree / strength
    s = sum(A,2);
    D = diag(s);

    % Handle zeros safely
    invsqrtD = diag( 1 ./ max(sqrt(s), eps) );

    % Normalized adjacency and Laplacian
    Ahat = invsqrtD * A * invsqrtD;       % eigenvalues in [-1, 1]
    Lsym = eye(size(A)) - Ahat;           % eigenvalues in [0, 2]
end
