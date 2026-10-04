function [yhat, eres] = reduced_fit_qr(y, Xb)
y = y(:);
[yhat,eres] = ols_projection([ones(numel(y),1),Xb],y);
end
