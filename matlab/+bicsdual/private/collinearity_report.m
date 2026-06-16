function colrep = collinearity_report(X, names, vif_cut, rho_cut)
if nargin<3, vif_cut = 10; end
if nargin<4, rho_cut = 0.98; end
names = ensure_names(names, size(X,2));
R = corr(X); R(1:size(R,1)+1:end)=0; % zero diag for reporting
[vif, bad_vif] = local_vif(X, vif_cut);
[ri, rj] = find(abs(R) > rho_cut);
pairs = arrayfun(@(k) sprintf('%s~%s (r=%.3f)', names{ri(k)}, names{rj(k)}, R(ri(k),rj(k))), 1:numel(ri), 'uni',0);
colrep = struct('vif',vif,'bad_vif_idx',find(bad_vif),'bad_vif_names',{names(bad_vif)}, ...
                'rho_gt_cut',{pairs},'rho_cut',rho_cut,'vif_cut',vif_cut);
end

function [vif, bad] = local_vif(X, cut)
% VIF with intercept
n = size(X,1);
Xc = [ones(n,1) X];
vif = nan(1,size(X,2));
for j=1:size(X,2)
    y = X(:,j);
    Z = Xc; Z(:,j+1)=[];  % regress y on remaining + intercept
    % QR solve
    beta = Z\y;
    yhat = Z*beta; res = y - yhat;
    R2 = 1 - (res'*res) / max(sum((y-mean(y)).^2), eps);
    vif(j) = 1 / max(1 - R2, 1e-12);
end
bad = vif > cut;
end

function names = ensure_names(names, p)
if nargin<1 || isempty(names), names = arrayfun(@(i) sprintf('X%02d',i),1:p,'uni',0); end
end
