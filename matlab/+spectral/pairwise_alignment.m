function result = pairwise_alignment(reference, target, kmax, do_sign_align)
%PAIRWISE_ALIGNMENT Direct caliber-versus-target spectral comparison.
%
% This reproduces legacy Module 1 while avoiding its manuscript-specific
% plotting dependencies. normalize_operators from topology_spectra is used.

if nargin < 4, do_sign_align = true; end
[Aref,~] = normalize_operators(reference);
[Atgt,~] = normalize_operators(target);
N = size(Aref,1);
kmax = min(kmax,N-1);

adj_ref = local_basis(Aref,kmax,'adj',false);
adj_tgt = local_basis(Atgt,kmax,'adj',false);

% Legacy Module 1 formed combinatorial Laplacians from normalized adjacency.
Lref = diag(sum(Aref,2)) - Aref;
Ltgt = diag(sum(Atgt,2)) - Atgt;
lap_ref = local_basis(Lref,kmax,'lap',true);
lap_tgt = local_basis(Ltgt,kmax,'lap',true);

result = struct();
result.adj = local_compare(adj_ref,adj_tgt,do_sign_align);
result.lap = local_compare(lap_ref,lap_tgt,do_sign_align);
end

function B = local_basis(M,k,mode,drop_first)
M = (M+M.')/2;
[V,e] = eig(full(M),'vector');
V = real(V);
e = real(e);
if strcmp(mode,'adj')
    [~,idx] = sort(abs(e),'descend');
else
    [~,idx] = sort(e,'ascend');
end
V = V(:,idx);
e = e(idx);
if drop_first
    V = V(:,2:k+1);
    e = e(2:k+1);
else
    V = V(:,1:k);
    e = e(1:k);
end
n = sqrt(sum(V.^2,1));
n(n<=eps)=1;
B.V = V ./ n;
B.eigenvalues = e;
end

function S = local_compare(R,X,do_sign_align)
K = min(size(R.V,2),size(X.V,2));
Vr = R.V(:,1:K);
Vx = X.V(:,1:K);
if do_sign_align
    s = sign(sum(Vr.*Vx,1));
    s(s==0)=1;
    Vx = Vx .* s;
end
C = min(abs(Vr.'*Vx),1);
diag_similarity = diag(C);
S = struct();
S.similarity_matrix = C;
S.diagonal_similarity = diag_similarity;
S.cumulative_similarity = cumsum(diag_similarity)./(1:K)';
S.eigenvalues_ref = R.eigenvalues(1:K);
S.eigenvalues_target = X.eigenvalues(1:K);
S.wasserstein = local_wasserstein(S.eigenvalues_ref,S.eigenvalues_target);
end

function w = local_wasserstein(a,b)
n = min(numel(a),numel(b));
w = mean(abs(sort(a(1:n))-sort(b(1:n))));
end
