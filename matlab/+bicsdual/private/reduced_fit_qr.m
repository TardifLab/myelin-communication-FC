function [yhat, eres] = reduced_fit_qr(y, Xb)
% Adds intercept, QR solve, returns fitted and residuals under reduced model.
y  = y(:); n = numel(y);
Xr = [ones(n,1) Xb];
[Qr, ~] = qr(Xr,0);
rb  = rank(Xr,1e-10); %#ok<NASGU>  % optional if you want
beta= Qr \ y;                       % economy solve (since Qr is orthonormal)
yhat= Xr * beta;
eres= y - yhat;
end
