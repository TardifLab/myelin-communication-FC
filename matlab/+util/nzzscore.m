function Z = nzzscore(X, dim)
%NZZSCORE Z-score nonzero finite entries while preserving matrix shape.
% dim=0 means vectorize finite, nonzero values; otherwise calls zscore-like
% behavior along the specified dimension while preserving NaNs/Infs as zeros.
if nargin < 2, dim = 0; end
X = double(X);
Z = X;
if dim == 0
    mask = isfinite(X) & X ~= 0;
    mu = mean(X(mask));
    sd = std(X(mask));
    if isempty(mu) || sd == 0 || ~isfinite(sd)
        Z(mask) = 0;
    else
        Z(mask) = (X(mask) - mu) ./ sd;
    end
    Z(~mask) = 0;
else
    X(~isfinite(X)) = NaN;
    mu = mean(X, dim, 'omitnan');
    sd = std(X, 0, dim, 'omitnan');
    Z = (X - mu) ./ sd;
    Z(~isfinite(Z)) = 0;
end
end
