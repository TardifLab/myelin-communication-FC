function [Dxfm, meta] = preprocess_network(W, cfg)
%PREPROCESS_NETWORK Reproduce the positive z-shift transformation used in the paper.
%
% Nonzero finite edges are z-scored and shifted to strictly positive values.
% Absent edges remain zero.

W = double(W);
N = size(W,1);
if size(W,2) ~= N
    error('Community input must be a square matrix.');
end

W(~isfinite(W)) = 0;
W(1:N+1:end) = 0;

asymmetry = max(abs(W - W.'), [], 'all');
if asymmetry > 1e-10
    switch lower(char(string(cfg.symmetry_mode)))
        case 'average'
            W = (W + W.') ./ 2;
        case 'upper'
            W = triu(W,1);
            W = W + W.';
        otherwise
            error('cfg.symmetry_mode must be ''average'' or ''upper''.');
    end
end

edge_mask = W ~= 0 & isfinite(W);
vals = W(edge_mask);
if isempty(vals)
    error('Community input contains no finite nonzero edges.');
end
sd = std(vals,0);
if ~isfinite(sd) || sd == 0
    error('Nonzero community edge weights have zero variance.');
end

Dxfm = zeros(size(W));
Dxfm(edge_mask) = (vals - mean(vals)) ./ sd;
min_nonzero = min(Dxfm(edge_mask));
shift = abs(min_nonzero);
if shift == 0
    shift = eps;
else
    shift = shift + 0.001 * shift;
end
Dxfm(edge_mask) = Dxfm(edge_mask) + shift;
Dxfm(1:N+1:end) = 0;

meta = struct();
meta.n_nodes = N;
meta.density = nnz(triu(edge_mask,1)) / (N*(N-1)/2);
meta.asymmetry_before = asymmetry;
meta.shift = shift;
meta.edge_mask = edge_mask;
end
