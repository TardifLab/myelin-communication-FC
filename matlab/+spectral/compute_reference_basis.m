function basis = compute_reference_basis(reference, kmax)
%COMPUTE_REFERENCE_BASIS Compute caliber adjacency and Laplacian bases.
%
% The normalized operators are built by the existing topology_spectra
% normalize_operators helper. Adjacency modes are ordered by descending
% absolute eigenvalue; Laplacian modes are ordered by ascending eigenvalue.

[Ahat, Lsym] = normalize_operators(reference);
N = size(reference,1);
kmax = min(kmax, N-1);

[Va, ea] = local_eigenpairs(Ahat, kmax, 'adj');
[Vl, el] = local_eigenpairs(Lsym, kmax, 'lap');

strength = sum(reference,2);
strength = max(strength, eps);
Dhalf = spdiags(strength.^0.5, 0, N, N);

basis = struct();
basis.Ahat = Ahat;
basis.Lsym = Lsym;
basis.Dhalf = Dhalf;
basis.adj.V = Va;
basis.adj.eigenvalues = ea;
basis.lap.V = Vl;
basis.lap.eigenvalues = el;
basis.Vbasis.adj.V = Va;
basis.Vbasis.lap.V = Vl;
end

function [V,e] = local_eigenpairs(M,k,mode)
M = (M + M.') ./ 2;
opts = struct('issym',true,'isreal',true);
try
    if strcmp(mode,'adj')
        [V,D] = eigs(M,k,'largestabs',opts);
    else
        [V,D] = eigs(M,k,'smallestabs',opts);
    end
    e = real(diag(D));
    V = real(V);
catch
    [V,D] = eig(full(M),'vector');
    e = real(D);
    V = real(V);
end

if strcmp(mode,'adj')
    [~,idx] = sort(abs(e),'descend');
else
    [~,idx] = sort(e,'ascend');
end
V = V(:,idx);
e = e(idx);
V = V(:,1:k);
e = e(1:k);
n = sqrt(sum(V.^2,1));
n(n<=eps) = 1;
V = V ./ n;
end
