% Exploratory analysis for paper 3
%
% Covers these portions of the analysis
%   1. Group-level data trends
%   2. Basic topology
%   3. Optimal community detection
%
% 2025 Mark C Nelson MNI
%--------------------------------------------------------------------------

%% Params
% Data
  MCstr         = 'MTsat-tm';                                               % {'MTsat' 'R1-tm' 'gratio-tm' 'gratio-ts'}
  MC2str         = 'gratio-ts';                                               % {'MTsat' 'R1-tm' 'gratio-tm' 'gratio-ts'}
  SCstr         = 'COMMITscl';                                              % {'COMMITscl' 'COMMITsclvol' 'COMMITsclNoS'};
  % SCstr         = 'COMMITsclNoS';                                              % {'COMMITscl' 'COMMITsclvol' 'COMMITsclNoS'};
  LoSstr        = 'COMMITsclLoS';                                           % filename for edge length
  parc          = 'schaefer-400';                                           % {'schaefer-200' 'schaefer-400'}
  thr           = 30;                                                       % CB threshold level
  str           = [' (' parc ', D=' num2str(thr) '%)'];                    % tacked onto plot labels
  strc          = [' (' parc ')'];                                         % tacked onto plot labels
% Options
  segstr        = 'cor';                                                    % connectome segmentation {'cor' 'sub' 'full'} 
  xfm           = 'log';                                                    % W-->L xfm {'inv' 'log'}
  bins          = 100;                                                      % histogram bins
  fntsz         = 20;                                                       % plot font sizes
  usefont       = 'Helvetica';
  txtopt        = {'Interpreter','tex'};                                    % supports special characters in titles
% Plot sizes
  pposcon2      = [-719 141 721 656];                                       % larger for connectomes with community labels
  pposhist      = [-2559 -123 561 420];                                      % histograms
  pposscat      = [-2560 -128 561 420];                                     % scatter plots
  pposcorr      = [-2559 193 1313 1015];                                    % correlations; alt: [-2559 867 2560 341];    
  pposviol      = [-2559 398 864 810];                                      % violin plots
  ppos_bar      = [-2559 878 404 330];
% Plot sizes for final figs
  pposcon_final = [-2559 788 483 420];

% Con & hist plot locations for MEG & BOLD data
  pposcon_all=[-2559 905 344 303; -2214 905 344 303;  -1869 905 344 303; -1524 905 344 303; -1179 905 344 303; -835 905 344 303; -492 905 344 303];
  pposhist_all=[-2559 522 344 303; -2214 522 344 303; -1869 522 344 303; -1524 522 344 303; -1179 522 344 303; -835 522 344 303; -492 522 344 303];
  pposcon2_all=[-2559 139 344 303; -2214 139 344 303; -1869 139 344 303; -1524 139 344 303; -1179 139 344 303; -835 139 344 303; -492 139 344 303];

% Plot locations for ddisp_mat outputs
  pposh1=[-2558 788 560 420; -1997  788 560 420; -1436 788 560 420; -875  788 560 420];
  pposh2=[-1443 530 721 656; -1443 -128 721 656; -721  530 721 656; -721 -128 721 656];
  pposh3=[-2558 367 560 420; -1997  367 560 420; -1436 367 560 420; -875  367 560 420];


% -------------------------------------------------------------------------
%% ================        SETUP       ====================== %
% -------------------------------------------------------------------------

% Paths
  runloc        = 'local'; 
  P             = a0_A3_pathDefs(runloc);                                 
  addpath(genpath(P.MWC))
  locLabel      = P.dataLabel;
  locDeriv      = P.dataDerivative;
% Useful objects
  cm            = colormap(smartcmaps('inferno'));
  cm2           = colormap(smartcmaps('bentcoolwarm')); 
  cm3           = colormap(smartcmaps('kindlmann')); 
  cm4           = colormap(slanCM('viridis'));                  close all;
  pinfo         = conn_getParcInfo({parc},locLabel,segstr);                 % structure with parcellation info
  pinfo2        = conn_getParcInfo({parc},locLabel,'cor');                  % for HCP data (no sctx)

% Main Colormaps
  CLRS=get_colors;
  cmap_fc=[0 0 0];                                                             % black
  cmap_ntwk=cmap_Yeo7Network;                                                  % Yeo 7-Networks
  % cmap_pred=[.90,.15,.15; .33,.53,.86; .80,.80,.80; .50,.50,.50];              % Predictors; COMMIT=red; MYLNPC1=blue; LoS=lightGray; ED=darkGray alts (red=.9047 .1918 .1988; blue=.2941 .5447 .7494)
  % cmap_pred=[CLRS.carrot; CLRS.aztcprpl; CLRS.livrplgreen; .8,.8,.8];            % Predictors
  cmap_pred=[CLRS.carrot; CLRS.aztcprpl; CLRS.livrplred; CLRS.livrplgreen; .8,.8,.8];            % Predictors
  cmap_5cm =[.55,.80,.99; .24,.74,.71; .64,.40,.84; .98,.51,.45; .95,.82,.18]; % 5-CM; SI=blue; SPE=green; NE=purple; DE=red; CMY=gold (close to Seguin2020)
  cmap_intract=[0.466 0.674 0.188; 0.9 0.6 0.6; 0.55 0.85 0.95];
  cmap_discon=[0 .447 .741; .85 .325 .098];                                    % DIS vs CONnected node pairs

% % % xx=size(tmp,1); myfig; for ii=1:xx; line([ii ii],ylim,'Color',tmp(ii,:),'LineWidth',15); hold on; end;

% % % % subjects
% % %   load([locDeriv '/3_subjects_sc-fc_mwc-final_' parc '.mat'], 'Ss');
% % %   rmsubs_names      = {'sub-03' 'sub-18' 'sub-18r' 'sub-20r' 'sub-21r' 'sub-22r' 'sub-23r' 'sub-24r' 'sub-25' 'sub-26' 'sub-27' 'sub-27r' 'sub-29r'};
% % %   HCex              = cellfun(@(c) find(strcmp(c,Ss)), rmsubs_names);
% % %   Ss(HCex)          = [];

% Metadata
  % Nsub      = length(Ss);
  Nstx      = pinfo.nsub;
  Nctx      = pinfo.ncor;
  Nnode     = Nctx + Nstx;
  Nntwk     = length(pinfo.clabels);
  maskut    = logical(triu(ones(Nnode)));                                   % mask for upper triangle WITH main diagonal
  maskmd    = logical(eye(Nnode));                                          % mask for main diagonal
  
% Shortform Network names
  LBL_ntwk=pinfo.clabels;
  LBL_ntwk{cellstrfind(LBL_ntwk,'Visual')}='VIS';
  LBL_ntwk{cellstrfind(LBL_ntwk,'SomMot')}='SMN';
  LBL_ntwk{cellstrfind(LBL_ntwk,'DorsAtt')}='DAN';
  LBL_ntwk{cellstrfind(LBL_ntwk,'SalVentAtt')}='VAN'; %    S/VAN    SN/VAN  
  LBL_ntwk{cellstrfind(LBL_ntwk,'Limbic')}='LIM';
  LBL_ntwk{cellstrfind(LBL_ntwk,'Control')}='CON';
  LBL_ntwk{cellstrfind(LBL_ntwk,'Default')}='DMN';


%% ================  Load & prep EDGE data 

  loaddir_edge           = [locDeriv '/0_finalData_edges/grp_' parc '_CBthr-' num2str(thr) '_distDep/'];
  
% --- Load EDGE data
  load([loaddir_edge SCstr '.mat'],'Dtg');                  SC=Dtg;         % Connection-Strength
  load([loaddir_edge MCstr '.mat'],'Dtg');                  MC=Dtg;         % Myelin
  load([loaddir_edge MC2str '.mat'],'Dtg');                 MC2=Dtg;        % g-ratio
  load([loaddir_edge LoSstr '.mat'],'Dtg');                 LoS=Dtg;        % Edge Length
  load([loaddir_edge 'FC.mat'],'Dtg');                      FC=Dtg;         % FC
  % load([locDeriv '/EuclideanDist_' parc '.mat'],'EuclDist');ED=EuclDist;    % Euclidean Distance
  ED=pinfo.eucliddist;
  clear Dtg

% --- Compute delays
  load([loaddir_edge 'gratio-ts.mat'],'Dtg');               GC=Dtg;         % g-ratio
  k=6.0;                                                                    % proportionality constant (Drakesmith et al.) (m/s per µm)
  d=2.5;                                                                    % axon diameter constant
  v=k*(d./GC); v(GC==0)=NaN;                                                % velocity matrix (m/s)
  % delays=(LoS./1000)./v; delays(GC==0)=0;                                   % delay matrix (s) : length (mm -> m) / velocity
  delays=LoS./v; delays(GC==0)=0;                                           % delay matrix (ms) : length (mm) / velocity
  r_delay=1./delays; r_delay(delays==0)=0;
  delaystr='delays';
  clear Dtg

% --- mask for non-zero edges
  masknz    = SC~=0;                                                       % mask for all non zero edges in connectomes (binary connectome NOT uniform across weights)

% ---- Normalize
  SCnrm     = SC      ./ ( max(SC(:))      + (max(SC(:))*.01)      );       % to avoid 0 lengths at existing edges
  MCnrm     = MC      ./ ( max(MC(:))      + (max(MC(:))*.01)      );
  MC2nrm    = MC2     ./ ( max(MC2(:))     + (max(MC2(:))*.01)     );
  delaysnrm = delays  ./ ( max(delays(:))  + (max(delays(:))*.01)  );
  EDnrm     = ED      ./ ( max(ED(:))      + (max(ED(:))*.01)      );
  LoSnrm    = LoS     ./ ( max(LoS(:))     + (max(LoS(:))*.01)     );
  % FCnrm     = FC   ./ ( max(FC(:))  + (max(FC(:))*.01) );


% Final network density
  d=sprintf('%1.f', density_und(SCnrm) * 100); tmp=strsplit(str,'=');  str=[tmp{1} '=' d '%)'];

% ---- Load HCP data (Shafiei 2022)
  loaddir2='~/Desktop/ConnectomeProj/myelinWeightedConnectome/2_modelingFC/data';
  megdir  = [loaddir2 '/0_compiledConns/megdata_shafiei22/'];
% % % % Schaefer 200
% % %   parcstr='schaefer200';
% % %   FCx_mri = readNPY([megdir 'groupFCmri_' parcstr '.npy']);
% % %   tmp=load([megdir 'groupFCmeg_aec_orth_' parcstr '.npy.mat']); FCx_meg_aec_orth=tmp.megfc; megbandlbls=tmp.bands;
% Schaefer 400
  parcstr='schaefer400'; Nnode2=400;
  tmp=load([megdir 'groupFCmeg_aec_orth_' parcstr '.npy.mat']);         FCx_meg=tmp.megfc; megbandstr=tmp.bands;
  FCx_mri = readNPY([megdir 'groupFCmri_' parcstr '.npy']);
  clear tmp
% Info
  xstr=[' (' parcstr ', Shafiei 2022)'];                                    % tacked onto plot labels
  Nbands=size(megbandstr,1);
  megbandlbl={}; for ii=1:Nbands; megbandlbl{ii}=strrep(megbandstr(ii,:),' ',''); end
% Mods: zero out main diagonal
  FCx_mri(1:Nnode2+1:end)=0;
  useD=FCx_meg;
  for ii=1:Nbands; D=squeeze(useD(ii,:,:)); D(1:Nnode2+1:end)=0; useD(ii,:,:)=D; end
  FCx_meg=useD;
% Filter mask for out of sample data
  masknz_oos=masknz(Nstx+1:end,Nstx+1:end);                                 % nz mask for out-of-sample data (accounts for stx)

% ---- Load Sensorimotor-Association Axis map (Sydnor 2021)
  saadir  = [loaddir2 '/0_compiledConns/liu23/'];
  SA_map  = readNPY([saadir 'archemap_axis_final_wh.npy']);


% % % % =================     Look at the data
% % %   opt={cm,pinfo,pposcon2,pposhist,bins,fntsz,0}; opt2={cm2,pinfo,pposcon2,pposhist,bins,fntsz,1};
% % %   D=SCnrm;     D(D==0)=nan; ddisp_mat(D,[SCstr str],     SCstr,opt{:});             % SC
% % %   D=MCnrm;     D(D==0)=nan; ddisp_mat(D,[MCstr str],     MCstr,opt{:});             % Myelin
% % %   D=delaysnrm; D(D==0)=nan; ddisp_mat(D,[delaystr str],  delaystr,opt{:});          % Delay
% % %   D=LoSnrm;    D(D==0)=nan; ddisp_mat(D,['LoS' str],     'LoS', opt{:});            % LoS
% % %   % D=EDnrm;  D(D==0)=nan; ddisp_mat(D,['ED, ' parc],'ED',  opt{:});                % ED
% % %   D=FC;        D(D==0)=nan; ddisp_mat(D,['FC' str],      'FC',  opt2{:});           % FC
% % % 
% % % 
% % % % -------- Plot external data (Shafiei 2022)
% % % % --- Histograms & Connectomes with sequential node ordering
% % %   ppostmp=getPlotPos([-2559 442 408 355],Nbands+1,0);
% % %   ppostmp2=getPlotPos([-2559 -58 408 420],Nbands+1,0);
% % % % BOLD rsfMRI
% % %   D=FCx_mri; tmpstr='rsfMRI'; D(isnan(D))=0; D(isinf(D))=0; Dnan=D; Dnan(Dnan==0)=nan; 
% % %   connplot(D,tmpstr,[],ppostmp(1,:),0,1,cm2); set(gca,'XTickLabel','','YTickLabel',''); font(fntsz,usefont);          
% % %   myfig([tmpstr xstr],ppostmp2(1,:)); histogram(Dnan,bins); title(tmpstr); font(fntsz,usefont);
% % % % MEG
% % %   useD=FCx_meg; tmpstr='MEG-AEC-orth'; 
% % %   for ii=1:Nbands; D=squeeze(useD(ii,:,:)); D(isnan(D))=0; D(isinf(D))=0; Dnan=D; Dnan(Dnan==0)=nan;  LBL=[tmpstr ' (' megbandlbl{ii} ')'];
% % %       connplot(D,LBL,[],ppostmp(ii,:),0,1,cm2); set(gca,'XTickLabel','','YTickLabel',''); font(fntsz,usefont);
% % %       myfig([LBL xstr],ppostmp2(ii,:)); histogram(Dnan,bins); xlabel(LBL); font(fntsz,usefont);
% % %   end
% % % 
% % % % --- Connectomes with nodes ordered by network
% % %   ppostmp=getPlotPos([-2559 237 618 560],Nbands,1);
% % % % BOLD rsfMRI
% % %   D=FCx_mri; tmpstr='rsfMRI'; D(isnan(D))=0; D(isinf(D))=0; Dnan=D; Dnan(Dnan==0)=nan; 
% % %   connplot(D,tmpstr,pinfo2,ppostmp(1,:),0,1,cm2);           
% % % % MEG: AEC Orth
% % %   useD=FCx_meg; tmpstr='MEG-AEC-orth';  
% % %   for ii=1:Nbands; D=squeeze(useD(ii,:,:)); D(isnan(D))=0; D(isinf(D))=0; Dnan=D; Dnan(Dnan==0)=nan;  LBL=[tmpstr ' (' megbandlbl{ii} ')'];
% % %       connplot(D,LBL,pinfo2,ppostmp(ii,:),0,1,cm2);           
% % %   end
% % % 
% % % 
% % % % =================  Plot data + Correlations (My Data)
% % % % Correlations: Edge weights with each other
% % %   ppostmp=getPlotPos([-2559 540 318 257],7,0); opt={str,cm3,ppostmp(1,:),fntsz};  %ppostmp=getPlotPos([-2559 527 360 270],7,0); 
% % % % --- Reference data
% % % % SC vs MC
% % %   X=SCnrm; Y=MCnrm; xlbl=[SCstr '-log10']; ylbl=MCstr; X(X==0)=nan; X=log10(X); X(isnan(X))=0;
% % %   ddisp_corr(X(Y~=0),Y(Y~=0),xlbl,ylbl,opt{:}); 
% % % 
% % % % Correlations: Edge Weights vs Edge Length
% % %   X=LoSnrm; xlbl='Edge Length'; opt={str,cm3,ppostmp(1,:),fntsz};
% % %   Y=SCnrm; ylbl=[SCstr '-log10']; Y(Y==0)=nan; Y=log10(Y); Y(isnan(Y))=0; ddisp_corr(X(Y~=0),Y(Y~=0),xlbl,ylbl,opt{:}); % LoS vs SC
% % %   Y=MCnrm; ylbl=MCstr; ddisp_corr(X(Y~=0),Y(Y~=0),xlbl,ylbl,opt{:}); % LoS vs Myelin
% % % 
% % % % Correlations: Edge Weights vs BOLD FC
% % %   Y=FC; ylbl='FC'; opt={str,cm3,ppostmp(1,:),fntsz};
% % %   X=SCnrm;   xlbl=[SCstr '-log10']; X(X==0)=nan; X=log10(X); X(isnan(X))=0; ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,opt{:}); % SC
% % %   X=MCnrm;   xlbl=MCstr; ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,opt{:});    % Myelin
% % % 
% % % % Correlatations: vs delays
% % %   Y=delaysnrm; ylbl=delaystr; opt={str,cm3,ppostmp(1,:),fntsz};
% % %   X=SCnrm;     xlbl=[SCstr '-log10']; X(X==0)=nan; X=log10(X); X(isnan(X))=0; ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,opt{:}); % SC
% % %   X=MCnrm;     xlbl=MCstr;         ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,opt{:});    % Myelin
% % %   X=LoSnrm;    xlbl='Edge Length'; ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,opt{:});    % length
% % %   X=FC;        xlbl='FC';          ddisp_corr(X(Y~=0),Y(Y~=0),xlbl,ylbl,opt{:});    % BOLD FC
% % % 
% % % % --- Correlations: delay & MC vs all FC
% % % % vs delay
% % %   X=delaysnrm; xlbl=delaystr; ppostmp=getPlotPos([-2559 527 360 270],7,0);
% % %   Y=FC; ylbl='FC'; ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,str,cm3,ppostmp(1,:),fntsz);    % BOLD FC
% % %   useD=FCx_meg;
% % %   for ii=1:Nbands; Y=squeeze(useD(ii,:,:)); Y(isnan(Y))=0; Y(isinf(Y))=0;  ylbl=megbandlbl{ii};
% % %       ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,str,cm3,ppostmp(ii+1,:),fntsz);    % MEG FC           
% % %   end
% % % % vs MC
% % %   X=MCnrm; xlbl=MCstr; ppostmp=getPlotPos([-2559 527 360 270],7,0);
% % %   Y=FC; ylbl='FC'; ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,str,cm3,ppostmp(1,:),fntsz);    % BOLD FC
% % %   useD=FCx_meg;
% % %   for ii=1:Nbands; Y=squeeze(useD(ii,:,:)); Y(isnan(Y))=0; Y(isinf(Y))=0;  ylbl=megbandlbl{ii};
% % %       ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,str,cm3,ppostmp(ii+1,:),fntsz);    % MEG FC           
% % %   end
% % % 
% % % 
% % % % % % % Correlations: RESIDUAL edge weights vs FC
% % % % % %   Cst={FC;  SCnrm; MCnrm; LoSnrm};
% % % % % %   Wst={'FC' SCstr  MCstr 'LoS'};
% % % % % %   [Crsd,Ws_rsd]=conn_controlW8(Cst,Wst,'los');
% % % % % %   opt={str,cm3,pposscat,fntsz};
% % % % % %   ii=1; Y=Crsd{ii}; ylbl= Ws_rsd{ii};
% % % % % %   ii=2; X=Crsd{ii}; xlbl=[Ws_rsd{ii} '-log10']; X(X==0)=nan; X=log10(X); X(isnan(X))=0; ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,opt{:}); % SC
% % % % % %   ii=3; X=Crsd{ii}; xlbl= Ws_rsd{ii}; ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,opt{:}); % Myelin


% -------------------------------------------------------------------------
  

%% --------   Group-level data trends (figs for paper)
  ppostmp=[-2558 568 283 229];
% --- vs Caliber
  dx=SC(masknz); dx=log10(dx); XLBL='caliber';
% vs MTsat
  dy=MCnrm(masknz); YLBL='MTsat'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont); 
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% vs g-ratio
  dy=MC2(masknz); YLBL='g-ratio'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% vs delay
  dy=delays(masknz); YLBL='delay (ms)'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% --- vs Edge Length
  dx=LoS(masknz); XLBL='edge length';
% vs MTsat
  dy=MCnrm(masknz); YLBL='MTsat'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% vs g-ratio
  dy=MC2(masknz); YLBL='g-ratio'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% vs delay
  dy=delays(masknz); YLBL='delay (ms)'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% --- vs BOLD FC
  dx=FC(masknz); XLBL='BOLD-FC';
% vs MTsat
  dy=MCnrm(masknz); YLBL='MTsat'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% vs g-ratio
  dy=MC2(masknz); YLBL='g-ratio'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% vs delay
  dy=delays(masknz); YLBL='delay (ms)'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% --- vs delay
  dx=delays(masknz); XLBL='delay (ms)';
% vs MTsat
  dy=MCnrm(masknz); YLBL='MTsat'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% vs g-ratio
  dy=MC2(masknz); YLBL='g-ratio'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
% --- vs g-ratio
  dx=MC2(masknz); XLBL='g-ratio';
% vs MTsat
  dy=MCnrm(masknz); YLBL='MTsat'; myfig('',ppostmp); heatscatter(dx,dy);
  ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
  title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));

% --- Including MySD
load([loaddir_edge 'MySD.mat'],'Dtg');                       mysd=Dtg;
dx=SC(masknz); dx=log10(dx); XLBL='caliber';
dy=mysd(masknz); dx(dy==0)=[]; dy(dy==0)=[];  YLBL='MySD'; myfig('',ppostmp); heatscatter(dx,dy);
ylabel(YLBL); xlabel(XLBL); colormap(cm3); axis square; font(fntsz-2,usefont);
title(sprintf('\\rho = %.2f', corr(dx,dy,'type','Spearman')));
%clim([2 300])



% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%% ==============        Weight-Length XFM         =============== %
% -------------------------------------------------------------------------

% Compute lengths
  switch xfm
    case 'log'
        if any(SCnrm(:)<0 | SCnrm(:)>1 | MCnrm(:)<0 | MCnrm(:)>1)
            error('connection-weights must be in the interval [0,1) to use the transform -log(w_ij) \n');
        end
        L_sc  = -log(SCnrm);    L_sc(isinf(L_sc))   = 0;
        L_mc  = -log(MCnrm);    L_mc(isinf(L_mc))   = 0;
        L_mc2 = -log(MC2nrm);   L_mc2(isinf(L_mc2)) = 0;
    case 'inv'
        L_sc  = 1./SCnrm;       L_sc(isinf(L_sc))   = 0;
        L_mc  = 1./MCnrm;       L_mc(isinf(L_mc))   = 0;
        L_mc2 = 1./MC2nrm;      L_mc2(isinf(L_mc2)) = 0;
  end
% Normalize lengths
  L_sc  = L_sc  ./ max(L_sc(:));
  L_mc  = L_mc  ./ max(L_mc(:));
  L_mc2 = L_mc2 ./ max(L_mc2(:));


% % % % ---- PLOT Lengths data
% % %   str1   = 'Lengths: ';
% % %   str2   = [', xfm=' xfm str];
% % %   str3   = 'L-';
% % % % Data Summary
% % %   opt={cm,pinfo,pposcon2,pposhist,bins,fntsz,0}; ppostmp=[-2559 307 1313 490];
% % %   D=L_sc; D(D==0)=nan; ddisp_mat(D,[str1 SCstr str2],[str3 SCstr],opt{:});  % SC
% % %   D=L_mc; D(D==0)=nan; ddisp_mat(D,[str1 MCstr str2],[str3 MCstr],opt{:});  % MC
% % % % Correlations: reference lengths networks
% % %   opt={str2,cm,pposscat,fntsz};
% % %   X=L_sc; Y=L_mc;  xlbl=[str3 SCstr]; ylbl=[str3 MCstr]; ddisp_corr(X(Y~=0),Y(Y~=0),xlbl,ylbl,opt{:}); % L-sc vs L-mc
% % % % Correlations: vs Edge Length
% % %   DPLT      = {LoSnrm L_sc L_mc}';
% % %   W8s       = {'Edge Length' [str3 SCstr] [str3 MCstr]};
% % %   pltstr    = ['Lengths-weighted networks' str];
% % %   plot_conn_vsw8(DPLT,W8s,'Edge Length',{parc},pltstr,'group'); set(gcf,'Position',ppostmp) % heat scatter vs LoS
% % % % Correlations: vs FC
% % %   DPLT      = {FC L_sc L_mc}';
% % %   W8s       = {'FC' [str3 SCstr] [str3 MCstr]};
% % %   pltstr    = ['Lengths-weighted networks' str];
% % %   plot_conn_vsw8(DPLT,W8s,'FC',{parc},pltstr,'group'); set(gcf,'Position',ppostmp) % heat scatter vs FC



% -------------------------------------------------------------------------
%% ==============       Basic Topology (EDGE data)       =============== %
% -------------------------------------------------------------------------

% ===============  Compute
% 1:5 degree, clustering, distance, GlobEffic_dist, GlobEffic_hops, 
% 6:7 mean weight, participation coeff
  Ntop=7;                                                                    
  BIN=double(SC~=0);                                                        % Binary network

% --- Binary Network
  D=BIN; tmpcell=cell(1,Ntop);
% Degree
  dt                = sum(D); 
  dt                = dt(:);      
  tmpcell{1}        = dt./max(dt(:));                                       % degree
% Clustering Coeff
  dt                = clustering_coef_bu(D);                                % clustering
  tmpcell{2}        = dt./max(dt(:));
% Distance & Hops
  [dt,dt2]          = distance_wei_floyd(D);                                % distance & hops 
  dt                = dt ./ max(dt(:));
  tmpcell{3}        = {dt dt2};
% Global efficiency (dist)
  [~,dt]            = charpath(dt);
  tmpcell{4}        = dt;
% Global efficiency (hops)
  [~,dt]            = charpath(dt2);
  tmpcell{5}        = dt;
% Mean Edge Weight 
% (N/A)
% Participation Coeff
  dt                = participation_coef(D,pinfo.cirois,0);
  tmpcell{7}        = dt;
% output
  top_bin           = tmpcell;

% --- SC
  D=SC; L=L_sc; tmpcell=cell(1,Ntop);
% Strength (weighted degree)
  dt                = sum(D); 
  dt                = dt(:);      
  tmpcell{1}        = dt./max(dt(:));
% Clustering Coeff
  dt                = clustering_coef_wu(D);
  tmpcell{2}        = dt./max(dt(:));
% Distance & Hops
  [dt,dt2]          = distance_wei_floyd(L);                                % distance & hops           
  dt(isinf(dt))     = 0;
  dt                = dt ./ max(dt(:));
  tmpcell{3}        = {dt dt2};
% Global efficiency (dist)
  [~,dt]            = charpath(dt);
  tmpcell{4}        = dt;
% Global efficiency (hops)
  [~,dt]            = charpath(dt2);
  tmpcell{5}        = dt;
% Mean Edge Weight
  dt                = tmpcell{1} ./ sum(BIN,2);
  dt                = dt ./ max(dt(:));
  tmpcell{6}        = dt;
% Participation Coeff
  dt                = participation_coef(D,pinfo.cirois,0);
  tmpcell{7}        = dt;
% output
  top_sc            = tmpcell;

% --- Myelin
  D=MC; L=L_mc; tmpcell=cell(1,Ntop);
% Strength (weighted degree)
  dt                = sum(D); 
  dt                = dt(:);      
  tmpcell{1}        = dt./max(dt(:));
% Clustering Coeff
  dt                = clustering_coef_wu(D);
  tmpcell{2}        = dt./max(dt(:));
% Distance & Hops
  [dt,dt2]          = distance_wei_floyd(L);                                % distance & hops           
  dt(isinf(dt))     = 0;
  dt                = dt ./ max(dt(:));
  tmpcell{3}        = {dt dt2};
% Global efficiency (dist)
  [~,dt]            = charpath(dt);
  tmpcell{4}        = dt;
% Global efficiency (hops)
  [~,dt]            = charpath(dt2);
  tmpcell{5}        = dt;
% Mean Edge Weight
  dt                = tmpcell{1} ./ sum(BIN,2);
  dt                = dt ./ max(dt(:));
  tmpcell{6}        = dt;
% Participation Coeff
  dt                = participation_coef(D,pinfo.cirois,0);
  tmpcell{7}        = dt;
% output
  top_mc            = tmpcell;

% --- Myelin 2 (g-ratio)
  D=MC2; L=L_mc2; tmpcell=cell(1,Ntop);
% Strength (weighted degree)
  dt                = sum(D); 
  dt                = dt(:);      
  tmpcell{1}        = dt./max(dt(:));
% Clustering Coeff
  dt                = clustering_coef_wu(D);
  tmpcell{2}        = dt./max(dt(:));
% Distance & Hops
  [dt,dt2]          = distance_wei_floyd(L);                                % distance & hops           
  dt(isinf(dt))     = 0;
  dt                = dt ./ max(dt(:));
  tmpcell{3}        = {dt dt2};
% Global efficiency (dist)
  [~,dt]            = charpath(dt);
  tmpcell{4}        = dt;
% Global efficiency (hops)
  [~,dt]            = charpath(dt2);
  tmpcell{5}        = dt;
% Mean Edge Weight
  dt                = tmpcell{1} ./ sum(BIN,2);
  dt                = dt ./ max(dt(:));
  tmpcell{6}        = dt;
% Participation Coeff
  dt                = participation_coef(D,pinfo.cirois,0);
  tmpcell{7}        = dt;
% output
  top_mc2           = tmpcell;

% --- Delays
  D=r_delay; L=delays; tmpcell=cell(1,Ntop);
% Strength (weighted degree)
  dt                = sum(D); 
  dt                = dt(:);      
  tmpcell{1}        = dt./max(dt(:));
% Clustering Coeff
  dt                = clustering_coef_wu(D);
  tmpcell{2}        = dt./max(dt(:));
% Distance & Hops
  [dt,dt2]          = distance_wei_floyd(L);                                % distance & hops           
  dt(isinf(dt))     = 0;
  dt                = dt ./ max(dt(:));
  tmpcell{3}        = {dt dt2};
% Global efficiency (dist)
  [~,dt]            = charpath(dt);
  tmpcell{4}        = dt;
% Global efficiency (hops)
  [~,dt]            = charpath(dt2);
  tmpcell{5}        = dt;
% Mean Edge Weight
  dt                = tmpcell{1} ./ sum(BIN,2);
  dt                = dt ./ max(dt(:));
  tmpcell{6}        = dt;
% Participation Coeff
  dt                = participation_coef(D,pinfo.cirois,0);
  tmpcell{7}        = dt;
% output
  top_delay            = tmpcell;


clear tmpcell
  

% ===============    Visualize

  D={top_bin top_sc top_mc top_mc2 top_delay}; TMPLBL1={'Binary' 'Caliber' 'MTsat' 'g-ratio' 'Delay/Rate'};
  ppostmp=getPlotPos([-2559 489 373 308],length(D),0);
% --- Surfaces 
% Degree
  ii=1; TMPLBL2='Degree'; 
  DPLT=cellfun(@(c) c{ii}, D,'UniformOutput',0);
  Hcx=plot_conn_surf(DPLT,pinfo2,'cortex',TMPLBL1,TMPLBL2); setsurf(Hcx,[0 1],cm4);
% Clustering
  ii=2; TMPLBL2='Clustering'; 
  DPLT=cellfun(@(c) c{ii}, D,'UniformOutput',0);
  Hcx=plot_conn_surf(DPLT,pinfo2,'cortex',TMPLBL1,TMPLBL2); setsurf(Hcx,[0 1],cm4);
% Mean Edge Weight
  ii=6; TMPLBL2='Mean Edge Weight'; 
  DPLT=cellfun(@(c) c{ii}, D(2:4),'UniformOutput',0);
  Hcx=plot_conn_surf(DPLT,pinfo2,'cortex',TMPLBL1(2:4),TMPLBL2); setsurf(Hcx,[0 1],cm4);
% Participation Coeff
  ii=7; TMPLBL2='Participation'; 
  DPLT=cellfun(@(c) c{ii}, D,'UniformOutput',0);
  Hcx=plot_conn_surf(DPLT,pinfo2,'cortex',TMPLBL1,TMPLBL2); setsurf(Hcx,[0 1],cm4);

% ---- SPL Connectomes
  ii=3;
  connplot(top_sc{ii}{1},   'Caliber SPL',pinfo); colormap(cm); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk); title ''; % caliber
  connplot(top_mc{ii}{1},   'MTsat SPL', pinfo);  colormap(cm); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk); title ''; % MTsat
  connplot(top_mc2{ii}{1},  'g-ratio SPL',pinfo); colormap(cm); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk); title ''; % g-ratio
  connplot(top_bin{ii}{1},  'Binary SPL', pinfo); colormap(cm); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk); title ''; % binary
  connplot(top_delay{ii}{1},'Delay/Rate SPL',  pinfo); colormap(cm); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk); title ''; % delay

% ----- Scatter plot: weighted vs binary
  ppostmp=getPlotPos([-2551 460 369 337],4,0);
% -- vs DEGREE
      ii=1; str3=' Degree (bin/wei)'; X=top_bin{ii}; xlbl='Binary';
    % Caliber
      Y=top_sc{ii}; ylbl='Caliber'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(1,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
        % title({sprintf('r = %.2f',r), ['\fontsize{18pt}\rm' sprintf('(p = %.3g)',pval)] }, 'Interpreter', 'tex');
        % text(min(X),max(Y),sprintf('r = %.2f, p = %.3g',r,pval),'FontSize',fntsz,'VerticalAlignment','top'); hold off;
    % Myelin
      Y=top_mc{ii}; ylbl='MTsat'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(2,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
    % Myelin 2
      Y=top_mc2{ii}; ylbl='g-ratio'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(3,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
    % Delay
      Y=top_delay{ii}; ylbl='Delay/Rate'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(4,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
% -- vs CLUSTERING
      ii=2; str3=' Clustering'; X=top_bin{ii}; xlbl='Binary';
    % Caliber
      Y=top_sc{ii}; ylbl='Caliber'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(1,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
    % Myelin
      Y=top_mc{ii}; ylbl='MTsat'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(2,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
    % Myelin 2
      Y=top_mc2{ii}; ylbl='g-ratio'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(3,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
    % Delay
      Y=top_delay{ii}; ylbl='Delay/Rate'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(4,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
% -- vs Participation
      ii=7; str3=' Participation'; X=top_bin{ii}; xlbl='Binary';
    % Caliber
      Y=top_sc{ii}; ylbl='Caliber'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(1,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
    % Myelin
      Y=top_mc{ii}; ylbl='Myelin'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(2,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 

% ----- Scatter plot: Myelin1 vs Caliber
  ylbl='MTsat'; ppostmp=getPlotPos([-2559 492 321 305],5,0); %ppostmp=getPlotPos([-2551 460 369 337],4,0);
% -- STRENGTH
      ii=1; str3=' Strength'; X=top_sc{ii}; xlbl='Caliber'; Y=top_mc{ii}; 
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(1,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal');
% -- MEAN EDGE WEIGHT
      ii=6; str3=' Mean Edge Weight'; X=top_sc{ii}; xlbl='Caliber'; Y=top_mc{ii}; 
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(2,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal');
% -- CLUSTERING
      ii=2; str3=' Clustering'; X=top_sc{ii}; xlbl='Caliber'; Y=top_mc{ii}; 
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(3,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
% -- DISTANCE
      ii=3; str3=' Distance'; X=top_sc{ii}{1}; Y=top_mc{ii}{1}; X=X(Y~=0); Y=Y(Y~=0); xlbl='Caliber'; ylbl='Myelin';
      p=polyfit(X,Y,1); c=polyval(p,X); [r,pval]=corr(X,Y);
      ddisp_corr(X,Y,xlbl,ylbl,[str3 str],cm3,ppostmp(4,:),fntsz);  font(fntsz,usefont);
      title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
% -- Participation
      ii=7; str3=' Participation'; X=top_sc{ii}; xlbl='Caliber'; Y=top_mc{ii};
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(5,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 

% ----- Scatter plot: Myelin2 vs Caliber
  ylbl='g-ratio'; ppostmp=getPlotPos([-2559 492 321 305],5,0); %ppostmp=getPlotPos([-2551 460 369 337],4,0);
% -- STRENGTH
      ii=1; str3=' Strength'; X=top_sc{ii}; xlbl='Caliber'; Y=top_mc2{ii}; 
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(1,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal');
% -- MEAN EDGE WEIGHT
      ii=6; str3=' Mean Edge Weight'; X=top_sc{ii}; xlbl='Caliber'; Y=top_mc2{ii}; 
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(2,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal');
% -- CLUSTERING
      ii=2; str3=' Clustering'; X=top_sc{ii}; xlbl='Caliber'; Y=top_mc2{ii}; 
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(3,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
% -- DISTANCE
      ii=3; str3=' Distance'; X=top_sc{ii}{1}; Y=top_mc2{ii}{1}; X=X(Y~=0); Y=Y(Y~=0); xlbl='Caliber'; ylbl='Myelin';
      p=polyfit(X,Y,1); c=polyval(p,X); [r,pval]=corr(X,Y);
      ddisp_corr(X,Y,xlbl,ylbl,[str3 str],cm3,ppostmp(4,:),fntsz);  font(fntsz,usefont);
      title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
% -- Participation
      ii=7; str3=' Participation'; X=top_sc{ii}; xlbl='Caliber'; Y=top_mc2{ii};
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(5,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 

% ----- Scatter plot: DELAY vs Caliber
  ppostmp=getPlotPos([-2559 492 321 305],5,0); %ppostmp=getPlotPos([-2551 460 369 337],4,0);
% -- STRENGTH
      ii=1; str3=' Strength'; X=top_sc{ii}; xlbl='Caliber'; Y=top_delay{ii}; ylbl='Delay';
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(1,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
% -- Mean Edge Weight
      ii=6; str3=' Mean Edge Weight'; X=top_sc{ii}; xlbl='Caliber'; Y=top_delay{ii}; ylbl='Delay';
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(2,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal');      
% -- CLUSTERING
      ii=2; str3=' Clustering'; X=top_sc{ii}; xlbl='Caliber'; Y=top_delay{ii}; ylbl='Delay';
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(3,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
% -- DISTANCE
      ii=3; str3=' Distance'; X=top_sc{ii}{1}; Y=top_delay{ii}{1}; X=X(Y~=0); Y=Y(Y~=0); xlbl='Caliber'; ylbl='Delay';
      p=polyfit(X,Y,1); c=polyval(p,X); [r,pval]=corr(X,Y);
      ddisp_corr(X,Y,xlbl,ylbl,[str3 str],cm3,ppostmp(4,:),fntsz);  font(fntsz,usefont);
      title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
% -- Participation
      ii=7; str3=' Participation'; X=top_sc{ii}; xlbl='Caliber'; Y=top_delay{ii}; ylbl='Delay';
      p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(5,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 


% --- Scatter plot: vs SA-axis
  ppostmp=getPlotPos([-2559 492 321 305],5,0); % ppostmp=getPlotPos([-2551 460 369 337],3,0);
% -- vs DEGREE
      ii=1; str3=' Degree (bin/wei)'; X=SA_map; xlbl='SA-axis';
    % Caliber
      Y=top_sc{ii}; ylbl='Caliber'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(1,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
    % Myelin1
      Y=top_mc{ii}; ylbl='MTsat'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(2,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
    % Myelin2
      Y=top_mc2{ii}; ylbl='g-ratio'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(3,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
    % Delay
      Y=top_delay{ii}; ylbl='Delay/Rate'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(4,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
    % Binary
      Y=top_bin{ii}; ylbl='Binary'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(5,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.3g',r,pval), ''},'FontWeight','normal'); 
% -- vs MEAN EDGE WEIGHT
      ii=6; str3=' Mean Edge Weight'; X=SA_map; xlbl='SA-axis';
    % Caliber
      Y=top_sc{ii}; ylbl='Caliber'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(1,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
    % Myelin1
      Y=top_mc{ii}; ylbl='MTsat'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(2,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
    % Myelin2
      Y=top_mc2{ii}; ylbl='g-ratio'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(3,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
    % Delay
      Y=top_delay{ii}; ylbl='Delay/Rate'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(4,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
% -- vs CLUSTERING
      ii=2; str3=' Clustering'; X=SA_map; xlbl='SA-axis';
    % Caliber
      Y=top_sc{ii}; ylbl='Caliber'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(1,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
    % Myelin1
      Y=top_mc{ii}; ylbl='MTsat'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(2,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
    % Myelin2
      Y=top_mc2{ii}; ylbl='g-ratio'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(3,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
    % Delay
      Y=top_delay{ii}; ylbl='Delay/Rate'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(4,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
    % Binary
      Y=top_bin{ii}; ylbl='Binary'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(5,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
% -- vs PARTICIPATION
      ii=7; str3=' Participation'; X=SA_map; xlbl='SA-axis';
    % Caliber
      Y=top_sc{ii}; ylbl='Caliber'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(1,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
    % Myelin
      Y=top_mc{ii}; ylbl='Myelin'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(2,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 
    % Binary
      Y=top_bin{ii}; ylbl='Binary'; p=polyfit(X,Y,1); yfit=polyval(p,X); [r,pval]=corr(X,Y);
      myfig([str3 str],ppostmp(3,:)); scatter(X,Y,'filled'); hold on; plot(X,yfit,'k-','LineWidth',2);
      xlabel(xlbl); ylabel(ylbl); font(fntsz,usefont); title({sprintf('r = %.2f, p = %.2e',r,pval), ''},'FontWeight','normal'); 



% --- Bar plots
  D={top_bin top_sc top_mc top_mc2 top_delay}; TMPLBL1={'Binary' 'Caliber' 'MTsat' 'g-ratio' 'Delay'};
  ppostmp1=[-2559 490 266 307]; ppostmp=getPlotPos(ppostmp1,5,1); ppostmp2=[-2547 490 243 307]; ppostmp3=[-1991 490 291 307]; 
  ppostmp4=[-2559 538 238 259]; ppostmp5=[-2559 199 192 259];
  cmaptmp=[.25 .25 .25; cmap_pred(1:4,:)]; Ntmp=length(D);
% mean hops
  ii=3; str3='Mean Hops'; xx=1;
  DPLT=cell2mat(cellfun(@(c) mean(mean(c{ii}{2})), D,'UniformOutput',0));
  myfig([str3 str],ppostmp(xx,:)); mod_B=bar(DPLT,0.6);  hold on;
  ylabel(str3); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
  mod_B.CData=cmaptmp; mod_B.FaceColor='flat'; set(gca,'XTickLabel',TMPLBL1,'XTickLabelRotation',45); 
% Global Efficiency (dist)
  ii=4; str3='Global Efficiency'; xx=2;
  DPLT=cell2mat(cellfun(@(c) c{ii}, D,'UniformOutput',0));
  myfig([str3 str],ppostmp(xx,:)); mod_B=bar(DPLT,0.6);  hold on;
  ylabel(str3); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
  mod_B.CData=cmaptmp; mod_B.FaceColor='flat'; set(gca,'XTickLabel',TMPLBL1,'XTickLabelRotation',45); 
% % % % Global Efficiency (hops)
% % %   ii=5; str3='GlobEffic (hops)'; xx=3;
% % %   DPLT=cell2mat(cellfun(@(c) c{ii}, D,'UniformOutput',0));
% % %   myfig([str3 str],ppostmp(xx,:)); mod_B=bar(DPLT,0.6);  hold on;
% % %   ylabel(str3); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %   mod_B.CData=cmaptmp; mod_B.FaceColor='flat'; set(gca,'XTickLabel',TMPLBL1,'XTickLabelRotation',45);
% --- MEAN EDGE WEIGHT vs SAaxis
  ii=6; str3='\rho'; str4='Mean Weight vs SAaxis'; X=SA_map; savecorr=zeros(1,Ntmp-1);
  for xx=2:Ntmp; Y=D{xx}{ii}; savecorr(xx-1)=corr(X,Y,'type','Pearson'); end
  myfig(str4,ppostmp3); mod_B=bar(savecorr,0.6);  hold on;
  ylabel(str3,'FontWeight','bold'); font(fntsz,usefont); grid on; set(gca,'XGrid','off'); ylim([-.99 .99]);
  mod_B.CData=cmaptmp(2:end,:); mod_B.FaceColor='flat'; set(gca,'XTickLabel',TMPLBL1(2:end),'XTickLabelRotation',45); 
  set(gca,'TickLength',[0 0]); % set(gcf,'Position',ppostmp4);
% --- CLUSTERING vs caliber
  ii=2; str3='\rho'; str4='Clustering vs caliber'; X=top_sc{ii}; savecorr=zeros(1,Ntmp-2);
  for xx=3:Ntmp; Y=D{xx}{ii}; savecorr(xx-2)=corr(X,Y,'type','Pearson'); end
  myfig(str4,ppostmp2); mod_B=bar(savecorr,0.6);  hold on;
  ylabel(str3,'FontWeight','bold'); font(fntsz,usefont); grid on; set(gca,'XGrid','off'); ylim([-.99 .99]);
  mod_B.CData=cmaptmp(3:end,:); mod_B.FaceColor='flat'; set(gca,'XTickLabel',TMPLBL1(3:end),'XTickLabelRotation',45); 
  set(gca,'TickLength',[0 0]); % set(gcf,'Position',ppostmp4); set(gcf,'Position',ppostmp5);
% --- MEAN EDGE WEIGHT vs caliber
  ii=6; str3='\rho'; str4='MeanEdgeWeight vs caliber'; X=top_sc{ii}; savecorr=zeros(1,Ntmp-2);
  for xx=3:Ntmp; Y=D{xx}{ii}; savecorr(xx-2)=corr(X,Y,'type','Pearson'); end
  myfig(str4,ppostmp2); mod_B=bar(savecorr,0.6);  hold on;
  ylabel(str3,'FontWeight','bold'); font(fntsz,usefont); grid on; set(gca,'XGrid','off'); ylim([-.99 .99]);
  mod_B.CData=cmaptmp(3:end,:); mod_B.FaceColor='flat'; set(gca,'XTickLabel',TMPLBL1(3:end),'XTickLabelRotation',45); 
  set(gca,'TickLength',[0 0]); % set(gcf,'Position',ppostmp4); set(gcf,'Position',ppostmp5);




% % % %% TEMP SECTION TO EVALUATE EDGE MYELIN WEIGHTS
% % % 
% % % w_raw = MC(:);   % vector of edge weights
% % % w_raw(w_raw==0)=nan;
% % % myfig; histogram(w_raw,100); title('Raw myelin weights');
% % % fprintf('min %.3g, median %.3g, max %.3g, IQR %.3g\n', nanmin(w_raw), nanmedian(w_raw), nanmax(w_raw), iqr(w_raw(~isnan(w_raw))));
% % % 
% % % 
% % % cv_global = nanstd(w_raw)/nanmean(w_raw);
% % % % per-node CV (assuming adjacency in sparse format: [i,j,w])
% % % % edges_i, edges_j, edges_w are vectors
% % % n = max(max(edges_i), max(edges_j));
% % % node_cv = zeros(n,1);
% % % for u=1:n
% % %     idx = (edges_i==u)|(edges_j==u);
% % %     vals = edges_w(idx);
% % %     if numel(vals)>1
% % %         node_cv(u) = std(vals)/mean(vals);
% % %     else
% % %         node_cv(u) = NaN;
% % %     end
% % % end
% % % fprintf('global CV = %.3g; median node CV = %.3g\n', cv_global, nanmedian(node_cv));
% % % 
% % % 
% % % 
% % % 
% % % %% Load your edge myelin adjacency matrix (400x400, symmetric, zeros on diagonal)
% % % % Example: myelin_adj = load('myelin_matrix.mat'); 
% % % % assume variable is myelin_adj
% % % 
% % %   useD=SC; Dstr='CALIBER';
% % % 
% % % % Ensure symmetry
% % %   useD = (useD + useD') ./ 2;
% % % 
% % % % Binary adjacency (presence/absence of edges)
% % %   binary_adj = useD > 0;
% % % 
% % % % Compute node degree and strength
% % %   degree_vec   = sum(binary_adj,2);           % binary degree
% % %   strength_vec = sum(useD,2);           % weighted strength
% % % 
% % % % Mean edge weight per node
% % %   mean_weight_vec = strength_vec ./ degree_vec;  % (NaN for isolated nodes)
% % % 
% % % % Summaries
% % %   fprintf('--- Edge %s Network Diagnostics ---\n',Dstr);
% % %   fprintf('Density = %.3f\n', nnz(binary_adj)/(numel(binary_adj)-size(binary_adj,1)));
% % %   fprintf('Min mean edge weight = %.4e\n', min(mean_weight_vec));
% % %   fprintf('Max mean edge weight = %.4e\n', max(mean_weight_vec));
% % %   fprintf('Std of mean edge weight = %.4e\n', std(mean_weight_vec));
% % % 
% % % % Correlation of strength vs. degree
% % %   [R,P] = corr(strength_vec, degree_vec);
% % %   fprintf('Correlation (strength vs. degree) = %.3f (p=%.3g)\n', R, P);
% % % 
% % % % Distribution diagnostics
% % %   myfig([Dstr ': ' str] ,[-2558 206 2287 591]); 
% % %   subplot(1,3,1); histogram(useD(useD>0),50);
% % %     xlabel(['Edge ' Dstr ' weight']); ylabel('Count'); title('Distribution of nonzero edge weights'); font(fntsz,usefont);
% % %   subplot(1,3,2); scatter(degree_vec,strength_vec,'filled');
% % %     xlabel('Degree'); ylabel('Strength'); title('Node strength vs degree'); font(fntsz,usefont);
% % %   subplot(1,3,3); histogram(mean_weight_vec,50);
% % %     xlabel('Mean edge weight per node'); ylabel('Count'); title('Distribution of mean edge weight'); font(fntsz,usefont);
% % % 
% % % 
% % % 
% % % %% Linear model myelin ~ caliber + length
% % %   SCxfm=mylog(SC,1);
% % %   D=[MC(:) SCxfm(:) LoS(:)]; D=D(all(D~=0,2),:);
% % %   y=D(:,1); 
% % %   X=D(:,2:3); X=[ones(size(X,1),1), X];
% % %   b = X \ y; yhat = X*b;
% % %   R2 = 1-sum((y-yhat).^2) / sum((y-mean(y)).^2);
% % %   % SSres=sum((y-yhat).^2); SStot=sum((y-mean(y)).^2); R2=1-SSres/SStot;
% % %   fprintf('R^2 of myelin ~ caliber + length = %.3f\n', R2);
% % % 
% % % % Assess how much of unexplained variance is noise
% % %   residuals = y - yhat;
% % %   var_res = var(residuals);                                % total residual variance
% % %   std_res = std(residuals);                                % standard deviation
% % % 
% % % % --- Residuals vs caliber & length
% % %   D=residuals;
% % % % Correlations
% % %   p=polyfit(D,X(:,2),1); yfit_1=polyval(p,D); [r_1,pval_1]=corr(D,X(:,2));
% % %   p=polyfit(D,X(:,3),1); yfit_2=polyval(p,D); [r_2,pval_2]=corr(D,X(:,3));
% % % % Plot
% % %   myfig(str,[-2558 308 2287 489]); 
% % %   subplot(1,3,1); histogram(D,50);
% % %     xlabel('Residuals'); ylabel('Count'); title('Residuals: myelin ~ 1 + caliber + length'); font(fntsz,usefont);
% % %   subplot(1,3,2); scatter(D,X(:,2),'filled'); hold on; plot(D,yfit_1,'k-','LineWidth',2);
% % %     xlabel('Residuals'); ylabel('Caliber'); title('Model residuals vs caliber'); font(fntsz,usefont);
% % %     text(min(D),max(X(:,2)),sprintf('r = %.2f, p = %.3g',r_1,pval_1),'FontSize',18,'VerticalAlignment','top'); hold off;
% % %   subplot(1,3,3); scatter(D,X(:,3),'filled'); hold on; plot(D,yfit_2,'k-','LineWidth',2);
% % %     xlabel('Residuals'); ylabel('Length'); title('Model residuals vs length'); font(fntsz,usefont);
% % %     text(min(D),max(X(:,3)),sprintf('r = %.2f, p = %.3g',r_2,pval_2),'FontSize',18,'VerticalAlignment','top'); hold off;
% % % 
% % % % --- Residuals vs binary degree & clustering
% % %   D=zeros(Nnode); D(masknz)=residuals; D=sum(D,2);
% % % % Correlations
% % %   p=polyfit(D,top_bin{1},1); yfit_deg=polyval(p,D); [r_deg,pval_deg]=corr(D,top_bin{1});
% % %   p=polyfit(D,top_bin{2},1); yfit_clu=polyval(p,D); [r_clu,pval_clu]=corr(D,top_bin{2});
% % % % Plot
% % %   myfig(str,[-2558 308 2287 489]); 
% % %   subplot(1,3,1); histogram(D,50);
% % %     xlabel('Sum of residuals'); ylabel('Count'); title('myelin ~ 1 + caliber + length'); font(fntsz,usefont);
% % %   subplot(1,3,2); scatter(D,top_bin{1},'filled'); hold on; plot(D,yfit_deg,'k-','LineWidth',2);
% % %     xlabel('Sum of residuals'); ylabel('Binary degree'); title('Model residuals vs degree'); font(fntsz,usefont);
% % %     text(min(D),max(top_bin{1}),sprintf('r = %.2f, p = %.3g',r_deg,pval_deg),'FontSize',18,'VerticalAlignment','top'); hold off;
% % %   subplot(1,3,3); scatter(D,top_bin{2},'filled'); hold on; plot(D,yfit_clu,'k-','LineWidth',2);
% % %     xlabel('Sum of residuals'); ylabel('Binary clustering'); title('Model residuals vs clustering'); font(fntsz,usefont);
% % %     text(min(D),max(top_bin{2}),sprintf('r = %.2f, p = %.3g',r_clu,pval_clu),'FontSize',18,'VerticalAlignment','top'); hold off;
% % % 
% % % 
% % % % --- Residuals vs FC
% % %   D=residuals; D2=FC(masknz); D3=FCx_meg(6,:,:); D3=D3(masknz);
% % % % Correlations
% % %   p=polyfit(D,D2,1); yfit_1=polyval(p,D); [r_1,pval_1]=corr(D,D2);
% % %   p=polyfit(D,D3,1); yfit_2=polyval(p,D); [r_2,pval_2]=corr(D,D3);
% % % % Plot
% % %   myfig(str,[-2558 308 2287 489]); 
% % %   subplot(1,3,1); histogram(D,50);
% % %     xlabel('Residuals'); ylabel('Count'); title('Residuals: myelin ~ 1 + caliber + length'); font(fntsz,usefont);
% % %   subplot(1,3,2); scatter(D,D2,'filled'); hold on; plot(D,yfit_1,'k-','LineWidth',2);
% % %     xlabel('Residuals'); ylabel('BOLD FC'); title('Model residuals vs BOLD FC'); font(fntsz,usefont);
% % %     text(min(D),max(D2),sprintf('r = %.2f, p = %.3g',r_1,pval_1),'FontSize',18,'VerticalAlignment','top'); hold off;
% % %   subplot(1,3,3); scatter(D,D3,'filled'); hold on; plot(D,yfit_2,'k-','LineWidth',2);
% % %     xlabel('Residuals'); ylabel('Gamma_h_i FC'); title('Model residuals vs Gamma_h_i FC'); font(fntsz,usefont);
% % %     text(min(D),max(D3),sprintf('r = %.2f, p = %.3g',r_2,pval_2),'FontSize',18,'VerticalAlignment','top'); hold off;
% % % 
% % % 
% % % 
% % % 
% % % %% Permutation test: do weights encode topology beyond degree? 
% % % % null: shuffle weights across existing edges and compute node strengths,
% % % % then compare correlation with degree for observed vs nulls
% % %   D=SC; Dstr='CALIBER'; nperm = 5000;
% % %   deg=sum(D>0,2); strength=sum(D,2);
% % %   [r_deg_str,p_deg_str]=corr(deg,strength,'Type','Spearman');
% % %   fprintf('Observed Spearman degree vs strength: r=%.3f (p=%.3g)\n', r_deg_str, p_deg_str);
% % %   r_perm=zeros(nperm,1);
% % %   edges=D(maskut==1 & masknz==1);
% % % % Permute & corr
% % %   for k = 1:nperm
% % %         disp(['Perm: ' num2str(k)]);
% % %         Dperm=zeros(Nnode);
% % %         Dperm(maskut==1 & masknz==1)=edges(randperm(numel(edges)));         % shuffle edge weights randomly across edges
% % %         Dperm = Dperm + Dperm';                                             % symmetrize
% % %         deg_perm=sum(Dperm>0, 2); str_perm=sum(Dperm, 2);                   % recompute degree & strength
% % %         r_perm(k)=corr(deg_perm,str_perm, 'Type','Spearman');               % null correlation
% % %   end
% % % % Empirical p-value
% % %   p_emp=mean(r_perm <= r_deg_str);
% % %   fprintf('Permutation test p=%.4f\n', p_emp);
% % % % Plot null distribution
% % %   myfig([Dstr ', nperm=' num2str(nperm) '; str']); histogram(r_perm,50); hold on
% % %   xline(r_deg_str,'r','LineWidth',3);
% % %   xlabel('Spearman \rho (deg vs strength)'); ylabel('Count');
% % %   title(sprintf('Permutation null (p=%.4f)',p_emp)); font(fntsz,usefont);



% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%% ============      Community Detection (EDGE data)       ============= %
% -------------------------------------------------------------------------

% --- Get log-xfm version of group SC
% Load subject data (better to log then group average?)
% % %   loaddir_edge_s = [locDeriv '/0_finalData_edges/sub_' parc '_CBthr-' num2str(thr) '_distDep/'];
% % %   load([loaddir_edge_s SCstr '.mat'],'Dts');
% % % % log-xfm 
% % %   D=mylog(Dts,1); D=groupavg(D,3); D=D./(max(D(:))+(max(D(:))*.01));  myfig; histogram(D(D~=0),100); title('log, group, norm');
% % %   D2=mylog(SC,1); D2(isnan(D2))=0; D2=D2./(max(D2(:))+(max(D2(:))*.01));  myfig; histogram(D2(D2~=0),100); title('group, log, norm');
% % %   D3=mylog(Dts,1); D3=nanmedian(D3,3); D3(isnan(D3))=0; D3=D3./(max(D3(:))+(max(D3(:))*.01));  myfig; histogram(D3(D3~=0),100); title('log, median-grp, norm');
  

% Settings
  ppostmp=getPlotPos([-2558 396 446 401],3,0);
  mod_Ngam      = 100;                                                      % sampling density for 1st pass
  mod_Ngam2     = 200;                                                      % sampling density for 2nd pass
  mod_Nreps     = 100;                                                      % Number of partitions for each gamma on pass 2
  GAM_BASE      = [0.1 5];                                                  % [0.1 10] for finer partitions

% --- Data
% MEG data
  mm=6; useD=squeeze(FCx_meg(mm,:,:)); pltstr=megbandlbl{mm}; modstr=[]; usecm=cm;
% Caliber
  % useD=SCnrm; pltstr=SCstr; modstr=[]; usecm=cm; CLIM=[prctile(useD(useD~=0),1) prctile(useD(useD~=0),99)];
% BOLD-FC
  % useD=FC; pltstr='BOLD-FC'; modstr='negative_sym'; usecm=cm2; CMAX=max(abs(useD(:))); CLIM=[CMAX*-1 CMAX];
% Myelin
  % useD=MCnrm; pltstr=MCstr; modstr=[]; usecm=cm; CLIM=[prctile(useD(useD~=0),1) prctile(useD(useD~=0),99)];
% Delay
  % useD=delays; pltstr=delaystr; modstr=[]; usecm=cm; CLIM=[prctile(useD(useD~=0),1) prctile(useD(useD~=0),99)];
% Rate
  % useD=r_delay; pltstr='Rate'; modstr=[]; usecm=cm; CLIM=[prctile(useD(useD~=0),1) prctile(useD(useD~=0),99)];
% Caliber (log-xfm)
  % D=mylog(SC,1); D(isnan(D))=0; D=D./(max(D(:))+(max(D(:))*.01)); useD=D; pltstr=[SCstr ' log']; usecm=cm; % log-SC
  % CLIM=[prctile(useD(useD~=0),1) prctile(useD(useD~=0),99)];
% prep data
D=useD; 
if ~isequal(D,D.'); D(tril(ones(Nnode))==1)=0; D=D+D.'; disp('SYMMETRY!'); end % enforce symmetry
D=nzzscore(D,0); MV=abs(min(D(:))); MV=MV+(.001*MV); D=D+MV;                   % z-score & shift to positive values
D(isnan(D))=0; Dxfm=D;                                                         % nans back to 0s

% ===== FIRST PASS   
% --- Define gamma range for 1st pass
% Data-informed approach
  % % % mod_gammin=min(Dxfm(Dxfm~=0)); mod_gammax=max(Dxfm(Dxfm~=0));
  % % % if abs(skewness(Dxfm(Dxfm~=0))) > 2
  % % %       mod_gamvals=exp(linspace(log(mod_gammin),log(mod_gammax), mod_Ngam)); % sample logarithmically
  % % %       % mod_gamvals=logspace(log(mod_gammin),log(mod_gammax),mod_Ngam);       % overshoots: 10^(ln(x))
  % % %       xopt={'Xscale','log'}; disp('log-scale');
  % % % else
  % % %       mod_gamvals=linspace(mod_gammin,mod_gammax,mod_Ngam);               % sample linearly
  % % %       xopt={'Xscale','linear'}; disp('linear-scale');
  % % % end
% log-spaced gamma over a fixed, data-agnostic range
  mod_gamvals = logspace(log10(GAM_BASE(1)), log10(GAM_BASE(2)), mod_Ngam);

  % parpool;
% --- Modularity: N=100 gamma over range min to max
  mod_ci=zeros(Nnode,mod_Ngam);                                             % community assignments
  for ii = 1:mod_Ngam

      % Using gamma as a resolution parameter
        mod_gam=mod_gamvals(ii); 
        [mod_ci(:,ii),~] = community_louvain(Dxfm,mod_gam,[],modstr);
        % [mod_ci(:,ii),~] = community_louvain(Dxfm,mod_gam,[],'negative_sym'); % FC < 0

      % Hard thresholding approach (destablized in skewed networks)
        % % % mod_B=Dxfm; mod_B(Dxfm<mod_gam)=0;
        % % % [mod_ci(:,ii),~] = community_louvain(mod_B,[],[],[]);             % all D > 0
        % % % % [mod_ci(:,ii),~] = community_louvain(mod_B,[],[],'negative_sym'); % FC < 0      
  end

% --- plot community assignments
  TMPSTR1=['Modularity: ' pltstr]; TMPSTR2='1st pass partition'; 
  opt={'XTick',[]}; 
  % opt={'XTick',1:mod_Ngam,'XTickLabelRotation',45};
% Plot 1: default clim
  myfig([TMPSTR1 str],ppostmp(1,:)); imagesc(mod_ci); xlabel('gammas'); ylabel('nodes'); title(TMPSTR2); 
  colorbar; colormap(cm); font(fntsz,usefont); set(gca,opt{:});
  xticklabels(arrayfun(@(x) sprintf('%.1g',x), mod_gamvals,'UniformOutput',0));
% Plot 2: mod clim
  myfig([TMPSTR1 str],ppostmp(2,:)); imagesc(mod_ci); xlabel('gammas'); ylabel('nodes'); title(TMPSTR2); 
  colorbar; colormap(cm); font(fntsz,usefont); set(gca,opt{:}); clim([0 10]);
  xticklabels(arrayfun(@(x) sprintf('%.1g',x), mod_gamvals,'UniformOutput',0));

% --- Get range for 2nd pass
% Plot number of CIs for each partition
  mod_ciu_count=max(mod_ci,[],1); 
  myfig([TMPSTR1 str],ppostmp(3,:)); plot(mod_gamvals,mod_ciu_count,'LineWidth',3); hold on;
  title('Range for 2nd pass'); font(fntsz,usefont); ylabel('N CIs'); xlabel('gammas');
% --- Find Yu & Yv
  Yu=find(mod_ciu_count>2,1,'first'); Yv=find(mod_ciu_count>(Nnode/2+1),1,'first');
% Double checks
  if isempty(Yu), Yu=1;end;     if isempty(Yv), Yv=numel(mod_gamvals);end   % Fallback if empty (avoid degenerate ranges)
  if Yu >= Yv
      % If the bracket collapses, widen to a safe default or just keep the coarse range
        Yu = max(1, floor(0.25*numel(mod_gamvals)));
        Yv = min(numel(mod_gamvals), ceil(0.90*numel(mod_gamvals)));
  end
  line([mod_gamvals(Yu) mod_gamvals(Yu)],ylim,'LineStyle','--','LineWidth',2,'Color','k');
  line([mod_gamvals(Yv) mod_gamvals(Yv)],ylim,'LineStyle','--','LineWidth',2,'Color','k');
  

% ===== SECOND PASS
% --- Define gamma range for 2nd pass
% Padding for added stability
  pad=2;
  mod_gammin=mod_gamvals(max(1,Yu-pad)); 
  mod_gammax=mod_gamvals(min(numel(mod_gamvals),Yv+pad));
  % mod_gammin=mod_gamvals(Yu); mod_gammax=mod_gamvals(Yv);                 % no padding
% % % % Original data-informed approach
% % %   if abs(skewness(Dxfm(Dxfm~=0))) > 2
% % %         mod_gamvals=exp(linspace(log(mod_gammin),log(mod_gammax), mod_Ngam2)); % sample logarithmically
% % %         xopt={'Xscale','log'}; disp('log-scale');
% % %   else
% % %         mod_gamvals=linspace(mod_gammin,mod_gammax,mod_Ngam2);              % sample linearly
% % %         xopt={'Xscale','linear'}; disp('linear-scale');
% % %   end
% log-spaced as in 1st level
  mod_gamvals=logspace(log10(mod_gammin), log10(mod_gammax), mod_Ngam2);

% Data storage
  mod_ci            = zeros(Nnode,mod_Ngam2,mod_Nreps);                      % community assignments
  mod_q             = zeros(mod_Ngam2,mod_Nreps);                            % modularity values
  mod_ciu           = zeros(Nnode,mod_Ngam2);                                % consensus communities

% --- Modularity: gamma2 in range gamma1(Yu) to gamma1(Yv)  
  for ii = 1:mod_Ngam2
        mod_gam = mod_gamvals(ii);
            
          % mod_B = Dxfm; mod_B(Dxfm<mod_gam) = 0;                            % hard thresholding approach
         
        fprintf('gam %1.2i ',ii);
        tic;
        parfor irep = 1:mod_Nreps
           % Resolution parameter approach
             [mod_ci(:,ii,irep),mod_q(ii,irep)] = community_louvain(Dxfm,mod_gam,[],modstr);
            % [mod_ci(:,ii,irep),mod_q(ii,irep)] = community_louvain(Dxfm,mod_gam,[],'negative_sym'); % FC < 0
          
           % % % % Hard thresholding approach
           % % %   [mod_ci(:,ii,irep),mod_q(ii,irep)] = community_louvain(mod_B,[],[],[]);             % all D > 0
           % % %   % [mod_ci(:,ii,irep),mod_q(ii,irep)] = community_louvain(mod_B,[],[],'negative_sym'); % FC < 0
        end
        
        mod_citemp      = squeeze(mod_ci(:,ii,:));
        mod_ag          = agreement(mod_citemp) / mod_Nreps;                    % agreement (probability)
        
        mod_cinull      = zeros(size(mod_citemp));
      % get null agreement
        parfor irep = 1:mod_Nreps
            mod_cinull(:,irep) = mod_citemp(randperm(Nnode),irep);
        end
        mod_agnull      = agreement(mod_cinull) / mod_Nreps;
        mod_tau         = mean(mod_agnull(:));
        %mod_tau = mean(mod_ag(:));
        
        % get consensus (see Lancichinetti & Fortunato (2012))
        mod_ciu(:,ii) = consensus_und(mod_ag,mod_tau,10);
        
        fprintf(' finished in %.2f s\n',toc);
  end


% --- plot consensus community assignments
  TMPSTR1=['Modularity: ' pltstr]; TMPSTR2='2nd pass consensus partitions'; 
  opt={'XTick',[]};  % opt=[{'XTick',[]} xopt{:}]; 
  % opt={'XTick',1:mod_Ngam,'XTickLabelRotation',45};
% Plot 1: default clim
  myfig([TMPSTR1 str],ppostmp(1,:)); imagesc(mod_ciu); xlabel('gammas'); ylabel('nodes'); title(TMPSTR2); 
  colorbar; colormap(cm); font(fntsz,usefont); set(gca,opt{:});
  xticklabels(arrayfun(@(x) sprintf('%.1g',x), mod_gamvals,'UniformOutput',0));
% Plot 2: mod clim
  myfig([TMPSTR1 str],ppostmp(2,:)); imagesc(mod_ciu); xlabel('gammas'); ylabel('nodes'); title(TMPSTR2); 
  colorbar; colormap(cm); font(fntsz,usefont); set(gca,opt{:}); clim([0 10]);
  xticklabels(arrayfun(@(x) sprintf('%.1g',x), mod_gamvals,'UniformOutput',0));
% Plot number of CIs for each partition
  mod_ciu_count=max(mod_ciu,[],1); 
  myfig([TMPSTR1 str],ppostmp(3,:)); plot(mod_gamvals,mod_ciu_count,'LineWidth',3); hold on;
  title('Consensus Partitions'); font(fntsz,usefont); ylabel('N CIs'); xlabel('gammas');

% --- Quick diagnostics
% How many isolates (strength ~ 0) appear? Should stay ~0 if you stopped thresholding
  Dt=sum(Dxfm,2); ppostmp2=getPlotPos([-2559 -539 446 401],2,0);
  myfig('',ppostmp2(1,:)); plot(mod_gamvals, arrayfun(@(ii) sum(Dt==0), 1:numel(mod_gamvals)),'LineWidth',3);
  set(gca,opt{:}); xlabel('\gamma'); ylabel('# isolates'); font(fntsz,usefont);
% Mean community size vs gamma
  myfig('',ppostmp2(2,:)); plot(mod_gamvals, Nnode ./ max(mean(mod_ci,3),[],1),'LineWidth',3);
  set(gca,opt{:}); xlabel('\gamma'); ylabel('Mean community size'); font(fntsz,usefont);


% ---------------- STABILITY SUMMARY (z-Rand mean/var) ----------------
% --- look at partition similarity
  mod_rzvals = zeros(mod_Ngam2,2);
  for ii = 1:mod_Ngam2
        mod_rz = zeros(mod_Nreps);
        mod_mask = triu(true(mod_Nreps), 1);
        for i = 1:(mod_Nreps - 1)
            parfor j = (i + 1):mod_Nreps
                % z-score of Rand index (similarity between partitions)
                 mod_rz(i,j) = fcn_randz(mod_ci(:,ii,i),mod_ci(:,ii,j));
            end
        end
        mod_rzvec = mod_rz(mod_mask);
        mod_rzvals(ii,1) = mean(mod_rzvec);
        mod_rzvals(ii,2) = var(mod_rzvec);
        fprintf('gam %i of %i\n',ii,mod_Ngam2);
  end

% --- Save Results
  M=[]; M.gamvals=mod_gamvals; M.ciu=mod_ciu; M.rzvals=mod_rzvals; M.useD=useD; M.Dxfm=Dxfm; M.Dstr=pltstr; M.q=mod_q; 
  svdir=[locDeriv '/x_modularity/grp_' parc];   
  if ~isfolder(svdir); mkdir(svdir); end
  save([svdir '/' pltstr '.mat'],'M')

% --- Load Results
  pltstr=SCstr;          usecm=cm;
  % pltstr='BOLD-FC';      usecm=cm2;
  % pltstr=MCstr;          usecm=cm;
  % pltstr='Rate';           usecm=cm;
  % pltstr=MC2str;          usecm=cm;
  % pltstr='delay';          usecm=cm;
  
  % mm=6; pltstr=megbandlbl{mm}; usecm=cm;
  
  
% Load
  svdir=[locDeriv '/x_modularity/grp_' parc];   
  load([svdir '/' pltstr '.mat'])
% Unpack
  mod_gamvals=M.gamvals; mod_ciu=M.ciu; mod_rzvals=M.rzvals; useD=M.useD; Dxfm=M.Dxfm; mod_q=M.q; 
% Settings
  TMPSTR1=['Modularity: ' pltstr]; TMPSTR2='2nd pass consensus partitions';
  CLIM=[prctile(useD(useD~=0),1) prctile(useD(useD~=0),99)];
% % % % Scale of D
% % %   if abs(skewness(Dxfm(Dxfm~=0))) > 2
% % %         xopt={'Xscale','log'}; disp('log-scale');
% % %   else
% % %         xopt={'Xscale','linear'}; disp('linear-scale');
% % %   end
% FC
  % CMAX=max(abs(useD(:))); CLIM=[CMAX*-1 CMAX];

% ==== Explore partitions

% - Select optimal partition algorithmically
  D=mod_rzvals; Dmean=D(:,1); Dvar=D(:,2); alpha=0.3;
  Dmean_norm        = (Dmean-min(Dmean))/(max(Dmean)-min(Dmean)); 
  Dvar_norm         = (Dvar-min(Dvar))/(max(Dvar)-min(Dvar));
  Dscore            = alpha*Dmean_norm + (1-alpha) * (1-Dvar_norm);
  [~,mod_ci_sel]    = max(Dscore)

% manually select partition with max mean & min var
  % % % mod_ci_sel=34; 
  tmp=[2.7833 2.7834]; mod_ci_sel=find(mod_gamvals>tmp(1) & mod_gamvals<tmp(2))
  % line([mod_gamvals(mod_ci_sel) mod_gamvals(mod_ci_sel)],ylim,'LineWidth', 3,'Color','k')


% --- set partition to use
  mod_cic   = mod_ciu(:,mod_ci_sel);
  Ncis_str  = ['Ncis = ' num2str(max(mod_cic))]
  mod_str   = sum(abs(Dxfm))';

% ------ Visualize
% --- Plot number of CIs for each partition
  mod_ciu_count=max(mod_ciu,[],1); thr=find(mod_ciu_count>50,1);
  myfig([TMPSTR1 str],[-2559 -539 446 401]); plot(mod_gamvals,mod_ciu_count,'LineWidth',3); hold on;
  title('Consensus Partitions'); font(fntsz,usefont); ylabel('N CIs'); xlabel('gammas');
% % % % draw lines for N CIs < 50
% % %   line([mod_gamvals(thr) mod_gamvals(thr)],[0 50],'LineWidth', 1,'Color','k','LineStyle','--');
% % %   line([0 mod_gamvals(thr)],[50 50],'LineWidth', 1,'Color','k','LineStyle','--');
% draw line for selected partition
  L=line([mod_gamvals(mod_ci_sel) mod_gamvals(mod_ci_sel)],ylim,'LineWidth', 3,'Color','k'); legend(L,['part# ' num2str(mod_ci_sel)])

% --- Line plot
  myfig([TMPSTR1 ': ' Ncis_str str],[-2559 396 1219 401]); 
  yyaxis left;  plot(mod_gamvals,mod_rzvals(:,1),'LineWidth',3); hold on; ylabel('partition mean');
  yyaxis right; plot(mod_gamvals,mod_rzvals(:,2),'LineWidth',3); hold on; ylabel('partition variance');
  xlabel('gammas'); font(fntsz,usefont); set(gca,'Xscale','log'); title(TMPSTR2); % legend({'mean randz','var randz'}); 
% % % % draw line for N CIs < 50
% % %   line([mod_gamvals(thr) mod_gamvals(thr)],ylim,'LineWidth', 1,'Color','k','LineStyle','--');
% draw line for selected partition
  line([mod_gamvals(mod_ci_sel) mod_gamvals(mod_ci_sel)],ylim,'LineWidth', 3,'Color','k','LineStyle','-');

% --- Visualize tradeoff
  myfig([TMPSTR1 str],[-1339 396 506 401]); scatter(Dvar, Dmean); 
  xlabel('variance'); ylabel('mean'); font(fntsz,usefont);
  text(Dvar(mod_ci_sel),Dmean(mod_ci_sel),['← partition:' num2str(mod_ci_sel)],'FontSize',18); 
  

% --- Plot matrix
  [~,i]     = sortrows([mod_cic mod_str],[1 -2]);
  myfig(str,[-832 396 488 401]); imagesc(useD(i,i)); colorbar; colormap(usecm); 
  title(Ncis_str); font(fntsz,usefont); axis square; clim(CLIM);

% % % % --- surface
% % %   TMPLBL1={['\bf \gamma = ' num2str(round(mod_gamvals(mod_ci_sel),2))]};
% % %   Hcx=plot_conn_surf(mod_cic,pinfo2,'cortex',TMPLBL1,TMPSTR1); 
% % %   setsurf(Hcx,[0 max(mod_cic)],cm4);

% --- 3D plot of communities
% 1st filter very small communities
  for ii = 1:max(mod_cic)
      if length(find(mod_cic==ii)) < 5
          mod_cic(mod_cic==ii)=nan;
      end
  end
% Renumber remaining
  remaining_cis=unique(mod_cic(~isnan(mod_cic))); Ntmp=length(remaining_cis);
  for ii=1:Ntmp; mod_cic(mod_cic==remaining_cis(ii))=ii; end
  disp(['N very small communities excluded: ' num2str(length(find(isnan(mod_cic))))])
% Then plot remainder
  coor=pinfo.coor; Nt=max(mod_cic); n=ceil(Nt^(1/3)); m=ceil(Nt/n);
  myfig([TMPSTR1 ';   partition:' num2str(mod_ci_sel)],[-2559 -275 1263 591]); 
  for ii = 1:max(mod_cic) 
        subplot(n,m,ii)
        scatter3(coor(:,1),coor(:,2),coor(:,3),50,[0.8 0.8 0.8]); hold on
        scatter3(coor(mod_cic==ii,1),coor(mod_cic==ii,2),coor(mod_cic==ii,3),50,'b','filled');
        view(90,90); axis image
  end
% --- surface
  TMPLBL1={['\bf \gamma = ' num2str(round(mod_gamvals(mod_ci_sel),2))]};
  Hcx=plot_conn_surf(mod_cic,pinfo2,'cortex',TMPLBL1,TMPSTR1); 
  setsurf(Hcx,[0 max(mod_cic)],cm4);


% ==================    Quantify stuff
  D=mod_cic;

% vs Yeo 7-rsNetworks
  Dt=[1:Nntwk]'; D2=nan(Nnode,1); for ii=1:Nntwk; D2(pinfo.cis{ii},1)=Dt(ii,1); end; YeoRSNs=D2;
  fcn_randz(D,YeoRSNs)

% --- Save/load all optimal CIs
% mod_cic_sc_8=mod_cic;
% mod_cic_sc_11=mod_cic;
% mod_cic_mc_6=mod_cic;
% mod_cic_mc_4=mod_cic;
% mod_cic_mc_7=mod_cic;
% mod_cic_mc2_4=mod_cic;
% mod_cic_mc2_8=mod_cic;
% mod_cic_mc2_9=mod_cic;
% mod_cic_delay_4=mod_cic;
% mod_cic_delay_5=mod_cic;
% mod_cic_delay_7=mod_cic;
% mod_cic_rate_7=mod_cic;
% mod_cic_fc_17=mod_cic;
% mod_cic_fc_4=mod_cic;
% mod_cic_fc_2=mod_cic;
% --- MEG
% mod_cic_delta_3=mod_cic;
% mod_cic_theta_3=mod_cic;
% mod_cic_alpha_3=mod_cic;
% mod_cic_beta_3=mod_cic;
% mod_cic_gammalo_3=mod_cic;
% mod_cic_gammalo_8=mod_cic;
% mod_cic_gammahi_3=mod_cic;
% mod_cic_gammahi_4=mod_cic;
% mod_cic_gammahi_5=mod_cic;
% mod_cic_gammahi_7=mod_cic;
% Save
  M=[]; 
  M.cic_sc_8=mod_cic_sc_8;  M.cic_sc_11=mod_cic_sc_11; 
  M.cic_mc_4=mod_cic_mc_4;   M.cic_mc_6=mod_cic_mc_6;  M.cic_mc_7=mod_cic_mc_7; 
  M.cic_mc2_4=mod_cic_mc2_4;   M.cic_mc2_8=mod_cic_mc2_8;  M.cic_mc2_9=mod_cic_mc2_9;
  M.cic_rate_7=mod_cic_rate_7;
  M.cic_delay_4=mod_cic_delay_4;   M.cic_delay_5=mod_cic_delay_5;   M.cic_delay_7=mod_cic_delay_7;
  M.cic_fc_2=mod_cic_fc_2;   M.cic_fc_4=mod_cic_fc_4;   M.cic_fc_17=mod_cic_fc_17;
  M.cic_delta_3=mod_cic_delta_3;
  M.cic_theta_3=mod_cic_theta_3;
  M.cic_alpha_3=mod_cic_alpha_3;
  M.cic_beta_3=mod_cic_beta_3;
  M.cic_gammalo_3=mod_cic_gammalo_3;  M.cic_gammalo_8=mod_cic_gammalo_8;
  M.cic_gammahi_3=mod_cic_gammahi_3;  M.cic_gammahi_4=mod_cic_gammahi_4;  M.cic_gammahi_5=mod_cic_gammahi_5; M.cic_gammahi_7=mod_cic_gammahi_7;
  M.lbls={'caliber-8' 'caliber-11'                  ... 
          'MTsat-4' 'MTsat-6'  'MTsat-7'            ... 
          'g-ratio-4' 'g-ratio-8'  'g-ratio-9'      ... 
          'rate-7'                                  ... 
          'delay-4' 'delay-5' 'delay-7'             ... 
          'BOLD-FC-2' 'BOLD-FC-4' 'BOLD-FC-17'      ...
          'delta-3'  'theta-3'  'alpha-3'  'beta-3' ...
          'gammalo-3'  'gammalo-8'                  ...
          'gammahi-3'  'gammahi-4'  'gammahi-5'  'gammahi-7'};
  svdir=[locDeriv '/x_modularity/grp_' parc];   
  if ~isfolder(svdir); mkdir(svdir); end
  save([svdir '/all_optimal_partitions.mat'],'M')
% Load
  svdir=[locDeriv '/x_modularity/grp_' parc];   
  load([svdir '/all_optimal_partitions.mat'])

% --- Bar plots
  D={M.cic_sc  M.cic_mc  M.cic_mc2  M.cic_fc_2  M.cic_fc_4}; Nd=length(D);
  TMPLBL1={'Caliber' 'MTsat' 'g-ratio' 'FC2' 'FC4'};
  ppostmp1=[-2559 490 358 307]; ppostmp2=[-2174 490 266 307];
  cmaptmp=[cmap_pred(1:3,:); .25 .25 .25; .5 .5 .5; ];
% Similarity with Yeo RSNs
  str3='RIz: Yeo RSNs'; xx=1; vals=nan(1,Nd);
  for ii=1:Nd; vals(ii)=fcn_randz(D{ii},YeoRSNs); end
  myfig([str3 str],ppostmp1); mod_B=bar(vals,0.6);  hold on;
  ylabel(str3); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
  mod_B.CData=cmaptmp; mod_B.FaceColor='flat'; set(gca,'XTickLabel',TMPLBL1,'XTickLabelRotation',45);
% Similarity with BOLD-FC CIs
  str3='RIz: BOLD-FC4'; xx=1; vals=nan(1,2);
  for ii=1:3; vals(ii)=fcn_randz(D{ii},M.cic_fc_4); end
  myfig([str3 str],ppostmp2); mod_B=bar(vals,0.6);  hold on;
  ylabel(str3); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
  mod_B.CData=cmaptmp(1:3,:); mod_B.FaceColor='flat'; set(gca,'XTickLabel',TMPLBL1(1:3),'XTickLabelRotation',45);




% -------------  Matrices of edge weights by CI

% --- Edge Length
  D=LoSnrm; usecm=cm; TMPSTR='Edge Length '; CLIM=[prctile(D(D~=0),1) prctile(D(D~=0),99)];
% caliber CIs
  mod_cic=M.cic_sc; TMPSTR2=[TMPSTR ', Caliber CIs'];
  mod_str=sum(abs(D))'; [~,i]=sortrows([mod_cic mod_str],[1 -2]);
  myfig(str,[-832 396 488 401]); imagesc(D(i,i)); colorbar; colormap(usecm); 
  title(TMPSTR2); font(fntsz,usefont); axis square; clim(CLIM);
% caliber CIs
  mod_cic=M.cic_mc; TMPSTR2=[TMPSTR ', Myelin CIs'];
  mod_str=sum(abs(D))'; [~,i]=sortrows([mod_cic mod_str],[1 -2]);
  myfig(str,[-832 396 488 401]); imagesc(D(i,i)); colorbar; colormap(usecm); 
  title(TMPSTR2); font(fntsz,usefont); axis square; clim(CLIM);
% BOLD CIs
  mod_cic=M.cic_fc_4; TMPSTR2=[TMPSTR ', FC CIs'];
  mod_str=sum(abs(D))'; [~,i]=sortrows([mod_cic mod_str],[1 -2]);
  myfig(str,[-832 396 488 401]); imagesc(D(i,i)); colorbar; colormap(usecm); 
  title(TMPSTR2); font(fntsz,usefont); axis square; clim(CLIM);

% --- BOLD-FC
  D=FC; usecm=cm2; TMPSTR='FC '; CMAX=max(abs(useD(:))); CLIM=[CMAX*-1 CMAX];
% caliber CIs
  mod_cic=M.cic_sc; TMPSTR2=[TMPSTR ', Caliber CIs'];
  mod_str=sum(abs(D))'; [~,i]=sortrows([mod_cic mod_str],[1 -2]);
  myfig(str,[-832 396 488 401]); imagesc(D(i,i)); colorbar; colormap(usecm); 
  title(TMPSTR2); font(fntsz,usefont); axis square; clim(CLIM);
% caliber CIs
  mod_cic=M.cic_mc; TMPSTR2=[TMPSTR ', Myelin CIs'];
  mod_str=sum(abs(D))'; [~,i]=sortrows([mod_cic mod_str],[1 -2]);
  myfig(str,[-832 396 488 401]); imagesc(D(i,i)); colorbar; colormap(usecm); 
  title(TMPSTR2); font(fntsz,usefont); axis square; clim(CLIM);
% BOLD CIs
  mod_cic=M.cic_fc_4; TMPSTR2=[TMPSTR ', FC CIs'];
  mod_str=sum(abs(D))'; [~,i]=sortrows([mod_cic mod_str],[1 -2]);
  myfig(str,[-832 396 488 401]); imagesc(D(i,i)); colorbar; colormap(usecm); 
  title(TMPSTR2); font(fntsz,usefont); axis square; clim(CLIM);



% -------------  Matrices of mean within/between edge weights by CI


% ------------- Boxplots of within edge weights for each CI


% ------------- Topology within CIs defined by caliber vs myelin?
% Clustering of FC?
% Strength of FC?
% Degree of SC?




%% ============      Community Detection (loop version)      ============= %
% SC also
% -------------------------------------------------------------------------  
  parpool;


% Settings
  mod_Ngam      = 100;                                                      % sampling density for 1st pass
  mod_Ngam2     = 200;                                                      % sampling density for 2nd pass
  mod_Nreps     = 100;                                                      % Number of partitions for each gamma on pass 2
  GAM_BASE      = [0.1 5];                                                  % [0.1 10] for finer partitions


for xxx=1:5

    switch xxx
        % --- Data
        % Caliber
          % useD=SCnrm; pltstr=SCstr; usecm=cm; CLIM=[prctile(useD(useD~=0),1) prctile(useD(useD~=0),99)];
        case 1
            % FC
              useD=FC; pltstr='BOLD-FC'; modstr='negative_sym';
        case 2
            % Myelin1
              useD=MCnrm; pltstr=MCstr; modstr=[];
        case 3
            % Myelin2
              useD=MC2nrm; pltstr=MC2str; modstr=[];
        case 4
            % Delay
              useD=delays; pltstr=delaystr; modstr=[];
        case 5
            % Rate
              useD=r_delay; pltstr='Rate'; modstr=[];
    end
    disp(['*-*-*-  RUNNING:   ' pltstr])

    % prep data
    D=useD; 
    if ~isequal(D,D.'); D(tril(ones(Nnode))==1)=0; D=D+D.'; disp('SYMMETRY!'); end % enforce symmetry
    D=nzzscore(D,0); MV=abs(min(D(:))); MV=MV+(.001*MV); D=D+MV;                   % z-score & shift to positive values
    D(isnan(D))=0; Dxfm=D;                                                         % nans back to 0s
    
    % ===== FIRST PASS   
    % --- Define gamma range for 1st pass
    % log-spaced gamma over a fixed, data-agnostic range
      mod_gamvals = logspace(log10(GAM_BASE(1)), log10(GAM_BASE(2)), mod_Ngam);
    
    % --- Modularity: N=100 gamma over range min to max
      mod_ci=zeros(Nnode,mod_Ngam);                                             % community assignments
      for ii = 1:mod_Ngam
    
          % Using gamma as a resolution parameter
            mod_gam=mod_gamvals(ii); 
            [mod_ci(:,ii),~] = community_louvain(Dxfm,mod_gam,[],modstr);

            % % % switch xxx
            % % %     case 1
            % % %         [mod_ci(:,ii),~] = community_louvain(Dxfm,mod_gam,[],'negative_sym'); % FC < 0 
            % % %     otherwise
            % % %         [mod_ci(:,ii),~] = community_louvain(Dxfm,mod_gam,[],[]);             % all D > 0
            % % % end
            %  
      end
        
    % --- Get range for 2nd pass
    % --- Find Yu & Yv
      mod_ciu_count=max(mod_ci,[],1); 
      Yu=find(mod_ciu_count>2,1,'first'); Yv=find(mod_ciu_count>(Nnode/2+1),1,'first');
    % Double checks
      if isempty(Yu), Yu=1;end;     if isempty(Yv), Yv=numel(mod_gamvals);end   % Fallback if empty (avoid degenerate ranges)
      if Yu >= Yv
          % If the bracket collapses, widen to a safe default or just keep the coarse range
            Yu = max(1, floor(0.25*numel(mod_gamvals)));
            Yv = min(numel(mod_gamvals), ceil(0.90*numel(mod_gamvals)));
      end
      
    
    % ===== SECOND PASS
    % --- Define gamma range for 2nd pass
    % Padding for added stability
      pad=2;
      mod_gammin=mod_gamvals(max(1,Yu-pad)); 
      mod_gammax=mod_gamvals(min(numel(mod_gamvals),Yv+pad));
    % log-spaced as in 1st level
      mod_gamvals=logspace(log10(mod_gammin), log10(mod_gammax), mod_Ngam2);
    % Data storage
      mod_ci            = zeros(Nnode,mod_Ngam2,mod_Nreps);                      % community assignments
      mod_q             = zeros(mod_Ngam2,mod_Nreps);                            % modularity values
      mod_ciu           = zeros(Nnode,mod_Ngam2);                                % consensus communities
    
    % --- Modularity: gamma2 in range gamma1(Yu) to gamma1(Yv)  
      for ii = 1:mod_Ngam2
            mod_gam = mod_gamvals(ii);
                
              % mod_B = Dxfm; mod_B(Dxfm<mod_gam) = 0;                            % hard thresholding approach
             
            fprintf('gam %1.2i ',ii);
            tic;
            parfor irep = 1:mod_Nreps
               % Resolution parameter approach
                 [mod_ci(:,ii,irep),mod_q(ii,irep)] = community_louvain(Dxfm,mod_gam,[],modstr);          
            end
            
            mod_citemp      = squeeze(mod_ci(:,ii,:));
            mod_ag          = agreement(mod_citemp) / mod_Nreps;                    % agreement (probability)
            
            mod_cinull      = zeros(size(mod_citemp));
          % get null agreement
            parfor irep = 1:mod_Nreps
                mod_cinull(:,irep) = mod_citemp(randperm(Nnode),irep);
            end
            mod_agnull      = agreement(mod_cinull) / mod_Nreps;
            mod_tau         = mean(mod_agnull(:));
            %mod_tau = mean(mod_ag(:));
            
            % get consensus (see Lancichinetti & Fortunato (2012))
            mod_ciu(:,ii) = consensus_und(mod_ag,mod_tau,10);
            
            fprintf(' finished in %.2f s\n',toc);
      end
       
    
    % ---------------- STABILITY SUMMARY (z-Rand mean/var) ----------------
    % --- look at partition similarity
      mod_rzvals = zeros(mod_Ngam2,2);
      for ii = 1:mod_Ngam2
            mod_rz = zeros(mod_Nreps);
            mod_mask = triu(true(mod_Nreps), 1);
            for i = 1:(mod_Nreps - 1)
                parfor j = (i + 1):mod_Nreps
                    % z-score of Rand index (similarity between partitions)
                     mod_rz(i,j) = fcn_randz(mod_ci(:,ii,i),mod_ci(:,ii,j));
                end
            end
            mod_rzvec = mod_rz(mod_mask);
            mod_rzvals(ii,1) = mean(mod_rzvec);
            mod_rzvals(ii,2) = var(mod_rzvec);
            fprintf('gam %i of %i\n',ii,mod_Ngam2);
      end
    
    % --- Save Results
      M=[]; M.gamvals=mod_gamvals; M.ciu=mod_ciu; M.rzvals=mod_rzvals; M.useD=useD; M.Dxfm=Dxfm; M.Dstr=pltstr; M.q=mod_q; 
      svdir=[locDeriv '/x_modularity/grp_' parc];   
      if ~isfolder(svdir); mkdir(svdir); end
      save([svdir '/' pltstr '.mat'],'M')  
end


%% ============      Community Detection (loop version)      ============= %
% MEG FC only
% -------------------------------------------------------------------------  
  parpool;


% Settings
  mod_Ngam      = 100;                                                      % sampling density for 1st pass
  mod_Ngam2     = 200;                                                      % sampling density for 2nd pass
  mod_Nreps     = 100;                                                      % Number of partitions for each gamma on pass 2
  GAM_BASE      = [0.1 5];                                                  % [0.1 10] for finer partitions
  modstr=[];

for xxx=5:6

    useD=squeeze(FCx_meg(xxx,:,:)); pltstr=megbandlbl{xxx};


    disp(['*-*-*-  RUNNING:   ' pltstr])

    % prep data
    D=useD; 
    if ~isequal(D,D.'); D(tril(ones(Nnode))==1)=0; D=D+D.'; disp('SYMMETRY!'); end % enforce symmetry
    D=nzzscore(D,0); MV=abs(min(D(:))); MV=MV+(.001*MV); D=D+MV;                   % z-score & shift to positive values
    D(isnan(D))=0; Dxfm=D;                                                         % nans back to 0s
    
    % ===== FIRST PASS   
    % --- Define gamma range for 1st pass
    % log-spaced gamma over a fixed, data-agnostic range
      mod_gamvals = logspace(log10(GAM_BASE(1)), log10(GAM_BASE(2)), mod_Ngam);
    
    % --- Modularity: N=100 gamma over range min to max
      mod_ci=zeros(Nnode,mod_Ngam);                                             % community assignments
      for ii = 1:mod_Ngam
    
          % Using gamma as a resolution parameter
            mod_gam=mod_gamvals(ii); 
            [mod_ci(:,ii),~] = community_louvain(Dxfm,mod_gam,[],modstr);

            % % % switch xxx
            % % %     case 1
            % % %         [mod_ci(:,ii),~] = community_louvain(Dxfm,mod_gam,[],'negative_sym'); % FC < 0 
            % % %     otherwise
            % % %         [mod_ci(:,ii),~] = community_louvain(Dxfm,mod_gam,[],[]);             % all D > 0
            % % % end
            %  
      end
        
    % --- Get range for 2nd pass
    % --- Find Yu & Yv
      mod_ciu_count=max(mod_ci,[],1); 
      Yu=find(mod_ciu_count>2,1,'first'); Yv=find(mod_ciu_count>(Nnode/2+1),1,'first');
    % Double checks
      if isempty(Yu), Yu=1;end;     if isempty(Yv), Yv=numel(mod_gamvals);end   % Fallback if empty (avoid degenerate ranges)
      if Yu >= Yv
          % If the bracket collapses, widen to a safe default or just keep the coarse range
            Yu = max(1, floor(0.25*numel(mod_gamvals)));
            Yv = min(numel(mod_gamvals), ceil(0.90*numel(mod_gamvals)));
      end
      
    
    % ===== SECOND PASS
    % --- Define gamma range for 2nd pass
    % Padding for added stability
      pad=2;
      mod_gammin=mod_gamvals(max(1,Yu-pad)); 
      mod_gammax=mod_gamvals(min(numel(mod_gamvals),Yv+pad));
    % log-spaced as in 1st level
      mod_gamvals=logspace(log10(mod_gammin), log10(mod_gammax), mod_Ngam2);
    % Data storage
      mod_ci            = zeros(Nnode,mod_Ngam2,mod_Nreps);                      % community assignments
      mod_q             = zeros(mod_Ngam2,mod_Nreps);                            % modularity values
      mod_ciu           = zeros(Nnode,mod_Ngam2);                                % consensus communities
    
    % --- Modularity: gamma2 in range gamma1(Yu) to gamma1(Yv)  
      for ii = 1:mod_Ngam2
            mod_gam = mod_gamvals(ii);
                
              % mod_B = Dxfm; mod_B(Dxfm<mod_gam) = 0;                            % hard thresholding approach
             
            fprintf('gam %1.2i ',ii);
            tic;
            parfor irep = 1:mod_Nreps
               % Resolution parameter approach
                 [mod_ci(:,ii,irep),mod_q(ii,irep)] = community_louvain(Dxfm,mod_gam,[],modstr);          
            end
            
            mod_citemp      = squeeze(mod_ci(:,ii,:));
            mod_ag          = agreement(mod_citemp) / mod_Nreps;                    % agreement (probability)
            
            mod_cinull      = zeros(size(mod_citemp));
          % get null agreement
            parfor irep = 1:mod_Nreps
                mod_cinull(:,irep) = mod_citemp(randperm(Nnode),irep);
            end
            mod_agnull      = agreement(mod_cinull) / mod_Nreps;
            mod_tau         = mean(mod_agnull(:));
            %mod_tau = mean(mod_ag(:));
            
            % get consensus (see Lancichinetti & Fortunato (2012))
            mod_ciu(:,ii) = consensus_und(mod_ag,mod_tau,10);
            
            fprintf(' finished in %.2f s\n',toc);
      end
       
    
    % ---------------- STABILITY SUMMARY (z-Rand mean/var) ----------------
    % --- look at partition similarity
      mod_rzvals = zeros(mod_Ngam2,2);
      for ii = 1:mod_Ngam2
            mod_rz = zeros(mod_Nreps);
            mod_mask = triu(true(mod_Nreps), 1);
            for i = 1:(mod_Nreps - 1)
                parfor j = (i + 1):mod_Nreps
                    % z-score of Rand index (similarity between partitions)
                     mod_rz(i,j) = fcn_randz(mod_ci(:,ii,i),mod_ci(:,ii,j));
                end
            end
            mod_rzvec = mod_rz(mod_mask);
            mod_rzvals(ii,1) = mean(mod_rzvec);
            mod_rzvals(ii,2) = var(mod_rzvec);
            fprintf('gam %i of %i\n',ii,mod_Ngam2);
      end
    
    % --- Save Results
      M=[]; M.gamvals=mod_gamvals; M.ciu=mod_ciu; M.rzvals=mod_rzvals; M.useD=useD; M.Dxfm=Dxfm; M.Dstr=pltstr; M.q=mod_q; 
      svdir=[locDeriv '/x_modularity/grp_' parc];   
      if ~isfolder(svdir); mkdir(svdir); end
      save([svdir '/' pltstr '.mat'],'M')  
end



%% Summarizing community structure

  % Yeo 7-Networks 
  Dt=[1:Nntwk]'; D2=nan(Nnode,1); for ii=1:Nntwk; D2(pinfo.cis{ii},1)=Dt(ii,1); end; YeoRSNs=D2;


% Load
  svdir=[locDeriv '/x_modularity/grp_' parc];   
  load([svdir '/all_optimal_partitions.mat'])

%% -------------------- CONFIG --------------------

% Optional: choose max number of contingency maps to show in Panel C
  MAX_CONTINGENCY = 7; ppostmp=[-2559 422 366 375]; ppostmp2=[-2559 1 1200 796];       % n < 10
  % MAX_CONTINGENCY = 12; ppostmp=[-2559 1 842 796]; ppostmp2=[-2559 1 1200 796];       % n < 15
  % MAX_CONTINGENCY = 20; ppostmp=[-2559 -408 1187 1205]; ppostmp2=[-2559 -392 1716 1189];        % n > 15
  

%% -------------------- COLLECT PARTITIONS --------------------
  fn = fieldnames(M);
  is_part = startsWith(fn, 'cic_');             % only partition fields
  part_fields = fn(is_part);

% If M.lbls exists and matches count, respect its ordering; otherwise use field order
  if isfield(M,'lbls') && numel(M.lbls) == numel(part_fields)
      % reorder part_fields to match M.lbls
        [ok, idx] = ismember(M.lbls, part_fields);
        if all(ok)
            part_fields = part_fields(idx);
            part_labels = M.lbls(:);
        else
            warning('Mismatch between M.lbls and cic_* fields; using field order.');
            part_labels = strrep(part_fields,'cic_','');
        end
  else
        part_labels = strrep(part_fields,'cic_',''); 
  end

% Manual tweaks to labels
  part_labels = strrep(part_labels,'_','-');
  part_labels = strrep(part_labels,'mc2','gratio');
  part_labels = strrep(part_labels,'mc','MTsat');
  part_labels = strrep(part_labels,'sc','caliber');
  part_labels = strrep(part_labels,'fc','BOLD');
  part_labels = strrep(part_labels,'gamma','gam');
  part_labels = strrep(part_labels,'hi','H');
  part_labels = strrep(part_labels,'lo','L');


% Concatenate partitions into C (N x P), ensure column vectors
  C = [];
  for i = 1:numel(part_fields)
        x = M.(part_fields{i});
        if size(x,2) > 1 && size(x,1) == 1, x = x.'; end
        C = [C, x(:)];
  end

% Remove some
  remove_data={'caliber-11' 'MTsat-6' 'MTsat-4' 'gratio-4' 'gratio-8' 'delay-4' ... 
               'delay-7' 'BOLD-2' 'BOLD-17' 'gamL-8' 'gamH-3' 'gamH-4' 'gamH-7'};
  remove_data=[remove_data 'gamL-3' 'gamH-5' 'delta' 'theta' 'alpha' 'beta' 'BOLD-4'];   % option to remove all FC data
  itmp=cell2mat(cellstrfind(part_labels,remove_data)); C(:,itmp)=[]; part_labels(itmp)=[];

% % % % Append Yeo-7 reference
% % %   C = [C, YeoRSNs(:)];
% % %   part_labels = [part_labels; {'Yeo-7'}];

% --- Basic checks
  N = size(C,1);
  if N ~= 400
        warning('Expecting 400 nodes; found N=%d. Proceeding anyway.', N);
  end
% Relabel each partition to 1..k compact, preserving order (stable)
  C = relabel_compact(C);
% Remove small communities (< 2)
  [C, filt_info] = filter_small_communities(C,2);
% Community counts per partition
  kvec = max(C,[],1).';
% Adjust labels to match filtered k
  for ii = 1:length(kvec)
        tmpstr=part_labels{ii};
        part_labels{ii}=[extractBefore(tmpstr,'-') '-' num2str(kvec(ii))];
  end

% Last checks
  Nd=length(part_labels);
  part_labels=cellfun(@(c) extractBefore(c,'-'), part_labels,'UniformOutput',0);



%% -------------------- PANEL B: AMI heatmap + VI dendrogram + k strip ----
% Compute pairwise similarities/distances using your function
  [AMI, VI] = part_similarity_matrix(C);

% Convert VI to condensed distance for clustering
% Guard against tiny numerical negatives on diagonal
  VI(1:size(VI,1)+1:end)=0;  D=squareform(VI,'tovector');
% Hierarchical clustering (average linkage) with optimal leaf order
  Z=linkage(D,'average');   ord=optimalleaforder(Z,D);
% Reorder matrices and labels
  AMI_ord=AMI(ord,ord);  VI_ord=VI(ord,ord);  labels_o=part_labels(ord);  kvec_o=kvec(ord);

% ----- Option to combine AMI & VI into single matrix
  A=AMI_ord;  V=VI_ord;
  V_sim=1-(V./max(V(:)));                                                   % Convert VI to similarity in [0,1], higher = more similar
% Build combined upper/lower matrix
  tM=nan(size(A));
  tM(triu(true(size(A)), 1)) = A(triu(true(size(A)), 1));                    % upper: AMI
  tM(tril(true(size(A)), -1)) = V_sim(tril(true(size(A)), -1));              % lower: 1 - VI/max
% Encode k (normalized) on main diag (or leave NaN)
  diag_mode = 'nan';  % 'nan' | 'k'
  switch diag_mode
        case 'k'
            k_norm = (kvec_o - min(kvec_o)) / max(1, (max(kvec_o)-min(kvec_o)));
            tM(1:size(tM,1)+1:end) = k_norm;                                  % diagonal in [0,1]
     case 'nan'
        % leave as NaN
  end

%% ---------- Plot: top dendrogram + aligned heatmap ----------
% Dendrogram on TOP

tmpfntsz=20; Tdim=3;
fig=myfig('Similarity of community partitions',ppostmp);

t=tiledlayout(fig,Tdim,Tdim,'TileSpacing','compact','Padding','compact');

% (Top) dendrogram
axD=nexttile(t,[1 Tdim]);
[Hd,~,~]=dendrogram(Z,0,'Orientation','top','Reorder',ord); axis off; font(tmpfntsz,usefont);
set(axD,'XTick',[],'YTick',[],'Box','off','Layer','top');  title(axD,'VI dendrogram');
% Thicken dendrogram lines
if ~iscell(Hd), Hd=num2cell(Hd); end; cellfun(@(h) set(h,'LineWidth',2.5),Hd);

% (Bottom) heatmap
axH=nexttile(t,[Tdim-1 Tdim]);  imagesc(axH,tM,[0 1]); % title(axH,'Partition similarity summary');
set(axH,'YDir','normal','Layer','top','Color','w'); font(tmpfntsz,usefont);
colormap(axH,cm);  cb=colorbar(axH);
cb.Label.String = 'Similarity (upper: AMI, lower: 1 - VI/max, diag: k norm)';
% set(axH,'XTick',1:Nd,'XTickLabel',labels_o,'XTickLabelRotation',45,'YTick',1:Nd,'YTickLabel',labels_o);
% set(axH,'XAxisLocation','top','XTick',1:Nd,'XTickLabel',labels_o,'XTickLabelRotation',60,'YTick',1:Nd,'YTickLabel',labels_o);
set(axH,'XAxisLocation','top','XTick',1:Nd,'XTickLabel',labels_o,'XTickLabelRotation',90,'YTick',1:Nd,'YTickLabel',labels_o);


% Hide NaNs (diagonal if you later switch it to NaN) with alpha
himg = findobj(axH,'Type','image');  set(himg,'AlphaData',~isnan(tM));

% Align axes widths 
drawnow;  % ensure positions are up-to-date before reading/writing them
posH = axH.Position;          % [left bottom width height]
posD = axD.Position;
posD(1) = posH(1);            % same left edge
posD(3) = posH(3);            % same width
set(axD,'Position',posD);

% Make their x-scales consistent (columns aligned with leaves)
  xlim(axH,[0.5 Nd+0.5]); xlim(axD,[0.5 Nd+0.5]); linkaxes([axD axH],'x');

% Optional: leaf labels on the dendrogram instead of heatmap, swap:
% set(axD,'XTick',1:n,'XTickLabel',labs,'XTickLabelRotation',45);
% set(axH,'XTickLabel',[]);

% Draw cell boxes (thin everywhere, thick on diagonal)
hold(axH,'on'); thinLW = 0.6;  diagLW = 2.6;
boxopt={'Clipping','on','HitTest','off','EdgeColor','k'};
for ii = 1:Nd
    for jj = 1:Nd
        rectangle(axH,'Position',[jj-0.5,ii-0.5,1,1],'LineWidth',thinLW,boxopt{:});
    end
end
for ii = 1:Nd; rectangle(axH,'Position',[ii-0.5,ii-0.5,1,1],'LineWidth',diagLW,boxopt{:}); end
hold(axH,'off');

% (Optional) overlay k values on the diagonal as text
showKtext = true;
if showKtext
    for ii = 1:Nd
        text(ii,ii,sprintf('%d',kvec_o(ii)),'HorizontalAlignment','center', ...
            'VerticalAlignment','middle','FontSize',tmpfntsz,'FontWeight','bold','Color','k');
    end
end

set(gcf,'Renderer','painters'); 
exportgraphics(gcf,'~/Downloads/3_VI-AMI-dendrogram-heatmap_nolbl.pdf','ContentType','vector');


%% -------------------- PANEL C: Align each partition to Yeo-7 ------------

% Append Yeo-7 reference
  C = [C, YeoRSNs(:)];
  part_labels = [part_labels; {'Yeo-7'}];

% Find reference 
  idx_yeo = cellstrfind(part_labels,'Yeo-7'); % find(strcmp(part_labels, 'Yeo-7'), 1, 'last');
  Yref    = C(:, idx_yeo);

% Test vectors
  P = size(C,2);
  show_list = setdiff(1:P, idx_yeo, 'stable');
  if numel(show_list) > MAX_CONTINGENCY
        show_list = show_list(1:MAX_CONTINGENCY);
  end

% Info
  % nplots = numel(show_list); ncols  = ceil(sqrt(nplots)); nrows  = ceil(nplots / ncols);
  nplots = numel(show_list); ncols  = nplots; nrows  = 1; ppostmp2=[-2173 480 1695 317];

% ----- Plot: Row-normalized (fraction or Yeo-r in each target community)
  myfig('Row-normalized contingency via Hungarian matching',ppostmp2);
  t2 = tiledlayout(nrows, ncols, 'TileSpacing','compact','Padding','compact');

for ii = 1:nplots
    j = show_list(ii);
    out = match_to_ref(Yref, C(:,j));     % struct with Cmap & col_order

    ax = nexttile(t2);

    % draw image
    imagesc(ax, out.Cmap_row(:, out.col_order), [0 1]);  axis image; 

    % tidy axes
    set(ax,'YDir','normal','Layer','top');
    R = size(out.Cmap_row,1);
    K = size(out.Cmap_row,2);

    % exact data limits for whole matrix
    xlim(ax,[0.5 K+0.5]);
    ylim(ax,[0.5 R+0.5]);

    % OPTIONAL: make cells visually closer to square without fighting tiles
    % (skip if you’re happy with default)
    % pbaspect(ax,[K R 1]);    % comment out if it causes vertical compression

    title(ax, part_labels{j});
    xlabel(ax,'Matched communities');
    colormap(ax, cm);
    set(ax,'XTick',[], ...
           'YTick',1:Nntwk, 'YTickLabel', LBL_ntwk);
    font(tmpfntsz, usefont);
end
title(t2,'Alignment to Yeo-7','FontSize',tmpfntsz+4,'FontName',usefont); 


% ----- Plot: Col-normalized (fraction or target-community coming from each Yeo-r)
  myfig('Col-normalized contingency via Hungarian matching',ppostmp2);
  t2 = tiledlayout(nrows, ncols, 'TileSpacing','compact','Padding','compact');

for ii = 1:nplots
    j = show_list(ii);
    out = match_to_ref(Yref, C(:,j));     % struct with Cmap & col_order

    ax = nexttile(t2);

    % draw image
    imagesc(ax, out.Cmap_col(:, out.col_order), [0 1]);  axis image; 

    % tidy axes
    set(ax,'YDir','normal','Layer','top');
    R = size(out.Cmap_col,1);
    K = size(out.Cmap_col,2);

    % exact data limits for whole matrix
    xlim(ax,[0.5 K+0.5]);
    ylim(ax,[0.5 R+0.5]);

    % OPTIONAL: make cells visually closer to square without fighting tiles
    % pbaspect(ax,[K R 1]);    % comment out if it causes vertical compression

    title(ax, part_labels{j});
    xlabel(ax,'Matched communities');
    colormap(ax, cm);
    set(ax,'XTick',[], ...
           'YTick',1:Nntwk, 'YTickLabel', LBL_ntwk);
    font(tmpfntsz, usefont);
end
title(t2,'Alignment to Yeo-7','FontSize',tmpfntsz+4,'FontName',usefont); 

set(gcf,'Renderer','painters'); 
exportgraphics(gcf,'~/Downloads/4_column-contingency-heatmap.pdf','ContentType','vector');



%% -------------------- PANEL C: Align each partition to caliber ------------

% Remove some
  itmp=cellstrfind(part_labels,'Yeo-7'); C(:,itmp)=[]; part_labels(itmp)=[];

% Find reference 
  idx_ref = cellstrfind(part_labels,'caliber');
  Yref    = C(:, idx_ref);
  LBL_caliber={'L-post' 'L-medl' 'L-latl' 'L-frnt' 'R-post' 'R-medl' 'R-latl' 'R-frnt'};
  Nlbl=length(LBL_caliber);

% Test vectors
  P = size(C,2);
  show_list = setdiff(1:P, idx_ref, 'stable');
  if numel(show_list) > MAX_CONTINGENCY
        show_list = show_list(1:MAX_CONTINGENCY);
  end

% Info
  % nplots = numel(show_list); ncols  = ceil(sqrt(nplots)); nrows  = ceil(nplots / ncols);
  nplots = numel(show_list); ncols  = nplots; nrows  = 1; ppostmp2=[-2173 480 1695 317];

% ----- Plot: Row-normalized (fraction of REF-r in each target community)
  myfig('Row-normalized contingency via Hungarian matching',ppostmp2);
  t2 = tiledlayout(nrows, ncols, 'TileSpacing','compact','Padding','compact');

for ii = 1:nplots
    j = show_list(ii);
    out = match_to_ref(Yref, C(:,j));     % struct with Cmap & col_order

    ax = nexttile(t2);

    % draw image
    imagesc(ax, out.Cmap_row(:, out.col_order), [0 1]);  axis image; 

    % tidy axes
    set(ax,'YDir','normal','Layer','top');
    R = size(out.Cmap_row,1);
    K = size(out.Cmap_row,2);

    % exact data limits for whole matrix
    xlim(ax,[0.5 K+0.5]);
    ylim(ax,[0.5 R+0.5]);

    % OPTIONAL: make cells visually closer to square without fighting tiles
    % (skip if you’re happy with default)
    % pbaspect(ax,[K R 1]);    % comment out if it causes vertical compression

    title(ax, part_labels{j});
    xlabel(ax,'Matched communities');
    colormap(ax, cm);
    set(ax,'XTick',[],'YTick',1:Nlbl,'YTickLabel',LBL_caliber);
    font(tmpfntsz, usefont);
end
title(t2,'Alignment to Caliber','FontSize',tmpfntsz+4,'FontName',usefont); 


% ----- Plot: Col-normalized (fraction or target-community coming from each REF-r)
  myfig('Col-normalized contingency via Hungarian matching',ppostmp2);
  t2 = tiledlayout(nrows, ncols, 'TileSpacing','compact','Padding','compact');

for ii = 1:nplots
    j = show_list(ii);
    out = match_to_ref(Yref, C(:,j));     % struct with Cmap & col_order

    ax = nexttile(t2);

    % draw image
    imagesc(ax, out.Cmap_col(:, out.col_order), [0 1]);  axis image; 

    % tidy axes
    set(ax,'YDir','normal','Layer','top');
    R = size(out.Cmap_col,1);
    K = size(out.Cmap_col,2);

    % exact data limits for whole matrix
    xlim(ax,[0.5 K+0.5]);
    ylim(ax,[0.5 R+0.5]);

    % OPTIONAL: make cells visually closer to square without fighting tiles
    % pbaspect(ax,[K R 1]);    % comment out if it causes vertical compression

    title(ax, part_labels{j});
    xlabel(ax,'Matched communities');
    colormap(ax, cm);
    set(ax,'XTick',[],'YTick',1:Nlbl,'YTickLabel',LBL_caliber);
    font(tmpfntsz, usefont);
end
title(t2,'Alignment to Caliber','FontSize',tmpfntsz+4,'FontName',usefont); 



%% -------------------- Surfaces --------------------

% -- SURFACES
TMPSTR1 = 'Community Partitions: '; cmaptmp=[.8 .8 .8; cm4];
for pp = 1 %1:Nd
    Hcx = plot_conn_surf(C(:,pp),pinfo2,'cortex',part_labels(pp),TMPSTR1);
    setsurf(Hcx,[0 kvec(pp)],cmaptmp);
end

set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/1-5-surf-gratio.pdf','ContentType','vector');

% COORDS
  coor=pinfo.coor;  TMPSTR1 = 'Community Partitions: '; 
  % Nmax=max(C(:)); n=ceil(Nmax^(1/3)); m=ceil(Nmax/n); ppostmp=[-2559 438 473 359];
  Nmax=max(C(:)); m=Nmax; n=1; ppostmp=[-2559 675 1514 122];
  for pp=1:5 
      DPLT=C(:,pp); Nt=max(DPLT); 
      myfig([TMPSTR1 ', ' part_labels{pp}],ppostmp); 
      for ii = 1:Nt 
            subplot(n,m,ii)
            scatter3(coor(:,1),coor(:,2),coor(:,3),50,[0.8 0.8 0.8]); hold on
            scatter3(coor(DPLT==ii,1),coor(DPLT==ii,2),coor(DPLT==ii,3),50,'b','filled');
            view(90,90); axis image; axis off;
      end
        set(gcf,'Renderer','painters'); 
        svstr=[num2str(pp) '-nodes-' part_labels{pp}];
        exportgraphics(gcf,['~/Downloads/2-' svstr '.pdf'],'ContentType','vector');
        close all; 
  end


%% ------------ SCRATCH

% % % %% -------------------- Plot: dendrogram + combined heatmap --------------------
% % % % Dendrogram on left
% % % tmpfntsz=18;
% % % myfig('Community detection summary',[-2559 1 1200 796]);
% % % t = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');
% % % 
% % % % (Left) VI dendrogram with labels
% % % nexttile(t,[4 1]);
% % % [~,~,perm] = dendrogram(Z, 0, 'Orientation','left', 'Reorder', 1:numel(ord));  %#ok<ASGLU>
% % % set(gca,'YDir','reverse','XTick',[],'Box','off'); title('VI dendrogram');
% % % yticks(1:Nd);   yticklabels(labels_o(perm));
% % % set(gca,'TickLabelInterpreter','none','FontSize',tmpfntsz); font(tmpfntsz,usefont);
% % % 
% % % 
% % % % (Right) Combined heatmap (upper=AMI, lower=1−VI/max, diag=k)
% % % nexttile(t,[4 3]); imagesc(tM,[0 1]);  axis image;  colormap(cm);  cb=colorbar;
% % % cb.Label.String = 'Similarity (upper: AMI, lower: 1 - VI/max, diag: k norm)';
% % % set(gca,'XTick',1:Nd,'XTickLabel',labels_o,'XTickLabelRotation',45,'YTick',1:Nd,'YTickLabel',labels_o,'Layer','top');
% % % title('Partition similarity summary'); font(tmpfntsz,usefont);
% % % 
% % % % Hide NaNs with alpha (if any)
% % % alphaMask = ~isnan(tM);
% % % himg = findobj(gca,'Type','image');
% % % set(himg, 'AlphaData', alphaMask);
% % % 
% % % % Draw black boxes (thin for all cells, thick on diagonal)
% % %   hold on; thinLW=0.5;  diagLW=2.5;     % box line widths
% % %   boxopt={'Clipping','on','HitTest','off','EdgeColor','k'};
% % %   for ii = 1:Nd, for jj = 1:Nd, rectangle('Position',[jj-0.5,ii-0.5,1,1],'LineWidth',thinLW,boxopt{:}); end, end
% % %   for ii = 1:Nd; rectangle('Position',[ii-0.5,ii-0.5,1,1],'LineWidth',diagLW,boxopt{:}); end
% % %   hold off;
% % % 
% % % % (Optional) overlay k values on the diagonal as text
% % % showKtext = true;
% % % if showKtext
% % %     for ii = 1:Nd
% % %         text(ii,ii,sprintf('%d',kvec_o(ii)),'HorizontalAlignment','center', ...
% % %             'VerticalAlignment','middle','FontSize',tmpfntsz,'FontWeight','bold','Color','k');
% % %     end
% % % end



