function [dR2,p] = safe_compare2(Xb,Xf,y)
% Can be used to switch between reg solvers: QR or fitlm
% Can't handle interactions in hierarchical design though!
%   - Old QR solver is not rank-revealing or nest-preserving
%   - Will return d<0
if isempty(y) || size(Xb,1) < (size(Xf,2)+5), dR2 = 0; p = 1; return; end
[dR2,p] = compare_models_fast(y, Xb, Xf);         % QR much faster
% [dR2,p] = compare_models(y, Xb, Xf);            % fitlm
if dR2 < 0 && dR2 > -1e-10, dR2 = 0; end          % numerical fluff
if ~isfinite(dR2), dR2 = 0; end
if ~isfinite(p), p = 1; end
end
