function [ALLCM,LBL,varargout] = communicationModeling_bin(BIN,DIST,varargin)
% Computes several communication models from BINARY structural 
% connectivity data.
% 
%       * Requires the Brain Connectivity Toolbox *
%
% Models include:
%   SPE (shortest paths efficiency) 
%   NE (navigation efficieny) 
%   SIE (inverse of search information)
%   PT (path transitivity)
%   CMY (communicability)
%   DE (diffusion efficiency; based on MFPT)
%
% Outputs incldue versions of these models which have been normalized 
% (z-scored & log-xfm if skewed) to support their use in regression
% modeling.
%
%
% Input
%                           + Required +
%   BIN     : NxN double; BINARY structural connectome
%   DIST    : NxN double; distance between nodes (e.g., Euclidean)
%                           + Optional +
%   SCstr   : char; label for SC to be appended to output CM labels
%
%
% Output
%                           + Required +
%   ALLCM   : 1x6 cell; all communication models (no transformations)
%   LBL     : 1x6 cell; labels
%                           + Optional +
%   ALLCMx  : 1x6 cell; transformed versions of communication models
%             (log-xfm if skewed and z-scored)
%   LBLx    : 1x6 cell; labels indicating transformations
%   
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------

%% Optional inputs
  nin      = max(nargin,1) - 2;
  defaults = {''};
  SCstr    = INhandler(varargin,nin,defaults);


%% Function

% Setup
  LBL={'SPE' 'NE' 'SIE' 'PT' 'CMY' 'DE'};                                   % labels for no xfm CMs

  if ~isempty(SCstr)
      LBL=cellfun(@(c) [SCstr '-' c],LBL,'UniformOutput',0);
  end

  LBLx=LBL;
  Ncm=length(LBL);
  ALLCMx=cell(1,Ncm); 
  ALLCM=cell(1,Ncm);
  Nnode=length(BIN);
  SKEWTOL=1.5;                                                              % limit of acceptable skewness      

% Create "lengths" version with infs for BIN=0
  L=BIN; L(L>0)=1; L(L==0)=inf;

% --- Shortest Path Efficiency
  xx=1; [D,~,~]=distance_wei_floyd(L);                                      % communication model
  D=1./D; D=(D+D.')/2; D(1:Nnode+1:end)=0; ALLCM{xx}=D;                     % inverse, symmetry average, null main diag and save
  if abs(skewness(D(~isinf(D)))) > SKEWTOL 
      D=mylog(D,0); LBLx{xx}=[LBLx{xx} 'log'];                              % log xfm
  end                   
  Dout=nzzscore(D,0); ALLCMx{xx}=Dout;                                      % zscore

% --- Navigation Efficiency
  xx=2; [~,D,~,~,~]=navigation_wu(BIN,DIST);                                % communication model
  D=1./D; D=max(D,D.'); D(1:Nnode+1:end)=0; ALLCM{xx}=D;                    % inverse, symmetry max, null main daig and save
  if abs(skewness(D(~isinf(D)))) > SKEWTOL 
      D=mylog(D,0); LBLx{xx}=[LBLx{xx} 'log'];                              % log xfm
  end
  Dout=nzzscore(D,0); ALLCMx{xx}=Dout;                                      % zscore

% --- Search Information
  xx=3; D=search_information(BIN,L);                                      % communication model
  D=1./D; D=(D+D.')/2; ALLCM{xx}=D;                                         % inverse, enforce symmetry and save (inverse non-standard for SI!)
  if abs(skewness(D(~isinf(D)))) > SKEWTOL 
      D=mylog(D,0); LBLx{xx}=[LBLx{xx} 'log'];                              % log xfm
  end
  Dout=nzzscore(D,0); ALLCMx{xx}=Dout;                                      % zscore

% --- Path Transitivity
  xx=4; D=path_transitivity(L); ALLCM{xx}=D;                              % communication model
  Dout=nzzscore(D,0); ALLCMx{xx}=Dout;                                      % zscore

% --- Communicability
  xx=5; S=diag(sum(BIN,2)); D=expm(S^(-1/2)*BIN*S^(-1/2));                  % communication model
  D=(D+D.')/2; D(1:Nnode+1:end)=0; ALLCM{xx}=D;                             % enforce symmetry, null main diagonal & save
  if abs(skewness(D(~isinf(D)))) > SKEWTOL 
      D=mylog(D,0); LBLx{xx}=[LBLx{xx} 'log'];                              % log xfm
  end
  Dout=nzzscore(D,0); ALLCMx{xx}=Dout;                                      % zscore

% --- Diffusion Efficiency
  xx=6; [~,D]=diffusion_efficiency(BIN);                                     % communication model
  D=(D+D.')/2; ALLCM{xx}=D;                                                 % enforce symmetry
  Dout=nzzscore(D,0); ALLCMx{xx}=Dout;                                      % zscore

LBLx=cellfun(@(c) [c 'z'],LBLx,'UniformOutput',0);

%% Optional output
nout = nargout-2;
if nout > 0
    varargout = cell(1,nout);
    for oo = 1 : nout
        switch oo
            case 1; varargout{oo}=ALLCMx;
            case 2; varargout{oo}=LBLx;
        end
    end
end

%--------------------------------------------------------------------------
end
