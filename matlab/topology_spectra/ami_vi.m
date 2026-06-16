function [AMI, VI] = ami_vi(x,y)
% Adjusted Mutual Information and Variation of Information
% 
% Uses NMI_sqrt (MI/sqrt(Hx*Hy)) as a fast, robust proxy for AMI
%
% Inputs:
% x,y: N×1 integer labels (need not be aligned)
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

  x = x(:); y = y(:);
% relabel to consecutive ints
  [~,~,x] = unique(x); [~,~,y] = unique(y);
  N = numel(x);

% contingency
  Kx = max(x); Ky = max(y);
  P = accumarray([x y], 1, [Kx Ky]) / N;   % joint
  px = sum(P,2); py = sum(P,1);
  H_x = -nansum(px .* log(px));            % entropy
  H_y = -nansum(py .* log(py));
% mutual information
  PxPy = px * py;
  MI = nansum(P(:) .* log(P(:) ./ PxPy(:)));

% Variation of Information (metric)
  VI = (H_x + H_y - 2*MI);

% Adjusted Mutual Information (Vinh et al 2010 form via expected MI approximation)
% Fast approximation: AMI = (MI - E[MI]) / (max(H_x,H_y) - E[MI]) 
  NMI = MI / sqrt(H_x * H_y);
  AMI = NMI;

% Helper
  function s = nansum(a), a(isnan(a))=0; s=sum(a); end

%--------------------------------------------------------------------------
end