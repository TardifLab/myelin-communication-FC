function score = compute_selection_score(R, alpha)
%COMPUTE_SELECTION_SCORE Combine high mean and low variance of z-Rand values.

if nargin < 2 || isempty(alpha), alpha = 0.30; end
mu = double(R.zrand_mean(:));
vr = double(R.zrand_variance(:));
valid = isfinite(mu) & isfinite(vr);
score = nan(size(mu));
if ~any(valid), return; end

mu_norm = normalize01(mu(valid));
vr_norm = normalize01(vr(valid));
score(valid) = alpha .* mu_norm + (1-alpha) .* (1-vr_norm);
end

function x = normalize01(x)
lo = min(x); hi = max(x);
if hi == lo
    x = zeros(size(x));
else
    x = (x-lo) ./ (hi-lo);
end
end
