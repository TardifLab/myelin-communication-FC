function [yhat, residual, r, Q] = ols_projection(X,y)
% SVD projection onto the entire column space, including redundant columns.
% Unit column norms keep rank decisions insensitive to predictor units.
X = full(double(X)); y = double(y(:));
scale = sqrt(sum(X.^2,1));
keep = scale>0;
X = X(:,keep); scale = scale(keep);
if isempty(X)
    Q=zeros(size(X,1),0); r=0; yhat=zeros(size(y)); residual=y; return;
end
X = bsxfun(@rdivide,X,scale);
[U,S,~] = svd(X,'econ'); s = diag(S);
tol = max(size(X))*eps(max(s));
r = sum(s>tol); Q=U(:,1:r);
yhat = Q*(Q'*y); residual=y-yhat;
end
