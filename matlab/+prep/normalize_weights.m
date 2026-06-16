function Xn = normalize_weights(X)
%NORMALIZE_WEIGHTS Normalize nonnegative edge weights to [0,~1).
X = double(X);
X(~isfinite(X)) = 0;
X(X < 0) = 0;
X = (X + X.') ./ 2;
X(1:size(X,1)+1:end) = 0;
mx = max(X(:));
if mx > 0
    Xn = X ./ (mx + mx*0.01);
else
    Xn = X;
end
end
