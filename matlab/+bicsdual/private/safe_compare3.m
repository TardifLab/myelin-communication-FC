function [dR2,p] = safe_compare3(Xb,Xf,y)
% Robust nested ΔR² + F-test.
% Includes intercept internally and handles rank-deficient designs using
% minimum-norm least-squares rather than rank-truncated raw columns.

dR2 = 0; 
p = 1;

if isempty(y) || isempty(Xb) || isempty(Xf)
    return
end

if size(Xb,1) ~= size(Xf,1) || size(Xb,1) ~= numel(y)
    return
end

y = y(:);
n = numel(y);

% Add intercept
Xb = [ones(n,1) Xb];
Xf = [ones(n,1) Xf];

% Drop non-finite rows defensively
ok = isfinite(y) & all(isfinite(Xb),2) & all(isfinite(Xf),2);
y  = y(ok);
Xb = Xb(ok,:);
Xf = Xf(ok,:);

n = numel(y);
if n < 3
    return
end

TSS = sum((y - mean(y)).^2);
if TSS <= eps
    return
end

% Numerical ranks
rb = rank(Xb);
rf = rank(Xf);

% If the full model adds no independent information, ΔR² is zero.
if rf <= rb
    dR2 = 0;
    p = 1;
    return
end

% Minimum-norm least-squares predictions.
% lsqminnorm handles rank-deficient designs without warning.
if exist('lsqminnorm','file') == 2
    betab = lsqminnorm(Xb, y);
    betaf = lsqminnorm(Xf, y);
else
    betab = pinv(Xb) * y;
    betaf = pinv(Xf) * y;
end

resb = y - Xb * betab;
resf = y - Xf * betaf;

rss_b = sum(resb.^2);
rss_f = sum(resf.^2);

% Enforce nesting monotonicity up to numerical precision
if rss_f > rss_b && rss_f - rss_b < 1e-10
    rss_f = rss_b;
end

dR2 = max(0, (rss_b - rss_f) / max(TSS, eps));

df1 = rf - rb;
df2 = n - rf;

if df1 <= 0 || df2 <= 0 || rss_f <= eps
    p = 1;
    return
end

F = ((rss_b - rss_f) / df1) / max(rss_f / df2, eps);

% Compute the upper-tail probability directly.
p = fcdf(F, df1, df2, 'upper');

% Avoid writing literal zeros to output tables due to numerical underflow.
if p == 0
    p = realmin;
end

if ~isfinite(p)
    p = 1;
end

end
