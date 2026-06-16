% Exploratory analysis for paper 3
%
% Covers these portions of the analysis
%   1. myelin communication vs caliber communication
%   2. communication efficiency distance relationship
%   3. spectral comparison of myelin vs caliber edge weights
%   4. Spectral fingerprinting of myelin and caliber communication
%
% 2025 Mark C Nelson MNI
%--------------------------------------------------------------------------

%% Params
% Data
  MCstr         = 'MTsat-tm';                                               % {'MTsat' 'R1-tm' 'gratio-tm' 'gratio-ts'}
  MC2str        = 'gratio-ts';                                               % {'MTsat' 'R1-tm' 'gratio-tm' 'gratio-ts'}
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
  cmap_5cm =[.24,.74,.71; .64,.40,.84; .55,.80,.99; .95,.82,.18; .98,.51,.45]; % 5-CM; SPE=green; NE=purple; SI=blue; CMY=gold DE=red; (close to Seguin2020)
  cmap_6cm =[.24,.74,.71; .64,.40,.84; .55,.80,.99; .59,.47,.45; .95,.82,.18; .98,.51,.45]; % 6-CM; SPE=green; NE=purple; SI=blue; PT=brown; CMY=gold DE=red; (close to Seguin2020)
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
  maskmd    = logical(eye(Nnode));                                          % mask for main diagonal maskmd=(repmat(1:Nnode,Nnode,1)==repmat((1:Nnode)',1,Nnode));
    

  
% Shortform Network names
  LBL_ntwk=pinfo.clabels;
  LBL_ntwk{cellstrfind(LBL_ntwk,'Visual')}='VIS';
  LBL_ntwk{cellstrfind(LBL_ntwk,'SomMot')}='SMN';
  LBL_ntwk{cellstrfind(LBL_ntwk,'DorsAtt')}='DAN';
  LBL_ntwk{cellstrfind(LBL_ntwk,'SalVentAtt')}='VAN'; %    S/VAN    SN/VAN  
  LBL_ntwk{cellstrfind(LBL_ntwk,'Limbic')}='LIM';
  LBL_ntwk{cellstrfind(LBL_ntwk,'Control')}='CON';
  LBL_ntwk{cellstrfind(LBL_ntwk,'Default')}='DMN';

% Yeo 7-Networks 
  Dt=[1:Nntwk]'; D2=nan(Nnode,1); for ii=1:Nntwk; D2(pinfo.cis{ii},1)=Dt(ii,1); end; YeoRSNs=D2;



%% ================  Load & prep EDGE data 

  loaddir_edge           = [locDeriv '/0_finalData_edges/grp_' parc '_CBthr-' num2str(thr) '_distDep/'];
  
% --- Load EDGE data
  load([loaddir_edge SCstr '.mat'],'Dtg');                  SC=Dtg;         % Connection-Strength
  load([loaddir_edge MCstr '.mat'],'Dtg');                  MC=Dtg;         % Myelin
  load([loaddir_edge MC2str '.mat'],'Dtg');                 MC2=Dtg;         % Myelin 2 (g-ratio)
  load([loaddir_edge LoSstr '.mat'],'Dtg');                 LoS=Dtg;        % Edge Length
  load([loaddir_edge 'FC.mat'],'Dtg');                      FC=Dtg;         % FC
  % load([locDeriv '/EuclideanDist_' parc '.mat'],'EuclDist');ED=EuclDist;    % Euclidean Distance
  ED=pinfo.eucliddist;

% --- Compute delays
  load([loaddir_edge 'gratio-ts.mat'],'Dtg');               GC=Dtg;         % g-ratio
  k=6.0;                                                                    % proportionality constant (Drakesmith et al.) (m/s per µm)
  d=2.5;                                                                    % axon diameter constant
  v=k*(d./GC); v(GC==0)=NaN;                                                % velocity matrix (m/s)
  delays=(LoS./1000)./v; delays(GC==0)=0;                                   % delay matrix (s) : length (mm -> m) / velocity
  r_delay=1./delays; r_delay(delays==0)=0;
  delaystr='Delay';
  clear Dtg
  
% ----- Load Community Partitions computed on SC & BOLD-FC data
  svdir=[P.dataDerivative '/x_modularity/grp_' parc];   
  load([svdir '/all_optimal_partitions.mat'])                   % M
  CIs_all=M;


% --- mask for non-zero edges
  masknz    = SC~=0;                                                       % mask for all non zero edges in connectomes (binary connectome NOT uniform across weights)

% ---- Normalize
  SCnrm     = SC   ./ ( max(SC(:))  + (max(SC(:))*.01)  );                  % to avoid 0 lengths at existing edges
  MCnrm     = MC   ./ ( max(MC(:))  + (max(MC(:))*.01)  );
  MC2nrm    = MC2  ./ ( max(MC2(:)) + (max(MC2(:))*.01)  );
  delaysnrm = delays  ./ ( max(delays(:))  + (max(delays(:))*.01)  );
  EDnrm     = ED   ./ ( max(ED(:))  + (max(ED(:))*.01)  );
  LoSnrm    = LoS  ./ ( max(LoS(:)) + (max(LoS(:))*.01) );
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


% ---- Combine all FC data (including z-scored version)
% Prep FC_bold IN sample: z-score
  D=FC; FCz_bold_in=nzzscore(D,0); FC_bold_in=D; 
  xx=1; tmplbls={'BOLDin'}; tmpyz={FCz_bold_in}; tmpy={FC_bold_in};
% Prep FC_bold OUT sample:z-score
  D=FCx_mri; FCz_bold_out=nzzscore(D,0); FC_bold_out=D; 
  xx=xx+1; tmplbls=[tmplbls 'BOLDout']; tmpyz=[tmpyz FCz_bold_out]; tmpy=[tmpy FC_bold_out];
% Prep FC_meg: z-score
  FCz_meg=cell(1,Nbands);  FC_meg=cell(1,Nbands);
  for ii = 1:Nbands
      D=squeeze(FCx_meg(ii,:,:)); FCz_meg{ii}=nzzscore(D,0); FC_meg{ii}=D;
  end
% Combine BOLD & MEG FC
  Nfc=Nbands+xx; ylbl=[tmplbls megbandlbl]; Yz_m=[tmpyz FCz_meg]; Y_m=[tmpy FC_meg];
  clear FCz_bold* FCz_meg FC_bold* FC_meg D dt Dt tmpmat D2 dt2

% Mod FC labels
  ylbl=strrep(ylbl,'BOLDin', 'BOLD_i_n');
  ylbl=strrep(ylbl,'BOLDout', 'BOLD_o_u_t');
  ylbl=strrep(ylbl,'lgamma', 'gamma_l_o');
  ylbl=strrep(ylbl,'hgamma', 'gamma_h_i');


% % % % =================     Look at the data
% % %   opt={cm,pinfo,pposcon2,pposhist,bins,fntsz,0}; opt2={cm2,pinfo,pposcon2,pposhist,bins,fntsz,1};
% % %   D=SCnrm;  D(D==0)=nan; ddisp_mat(D,[SCstr str],   SCstr,opt{:});             % SC
% % %   D=MCnrm;  D(D==0)=nan; ddisp_mat(D,[MCstr str],   MCstr,opt{:});             % Myelin
% % %   D=LoSnrm; D(D==0)=nan; ddisp_mat(D,['LoS' str],  'LoS', opt{:});             % LoS
% % %   % D=EDnrm;  D(D==0)=nan; ddisp_mat(D,['ED, ' parc],'ED',  opt{:});             % ED
% % %   D=FC;     D(D==0)=nan; ddisp_mat(D,['FC' str],   'FC',  opt2{:});            % FC
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
% % %   opt={str,cm3,pposscat,fntsz};
% % % % --- Reference data
% % % % SC vs MC
% % %   X=SCnrm; Y=MCnrm; xlbl=[SCstr '-log10']; ylbl=MCstr; X(X==0)=nan; X=log10(X); X(isnan(X))=0;
% % %   ddisp_corr(X(Y~=0),Y(Y~=0),xlbl,ylbl,opt{:}); 
% % % 
% % % % Correlations: Edge Weights vs Edge Length
% % %   X=LoSnrm; xlbl='Edge Length'; opt={str,cm3,pposscat,fntsz};
% % %   Y=SCnrm; ylbl=[SCstr '-log10']; Y(Y==0)=nan; Y=log10(Y); Y(isnan(Y))=0; ddisp_corr(X(Y~=0),Y(Y~=0),xlbl,ylbl,opt{:}); % LoS vs SC
% % %   Y=MCnrm; ylbl=MCstr; ddisp_corr(X(Y~=0),Y(Y~=0),xlbl,ylbl,opt{:}); % LoS vs Myelin
% % % 
% % % % Correlations: Edge Weights vs FC
% % %   Y=FC; ylbl='FC'; opt={str,cm3,pposscat,fntsz};
% % %   X=SCnrm; xlbl=[SCstr '-log10']; X(X==0)=nan; X=log10(X); X(isnan(X))=0; ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,opt{:}); % SC
% % %   X=MCnrm; xlbl=MCstr; ddisp_corr(X(X~=0),Y(X~=0),xlbl,ylbl,opt{:}); % Myelin
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
  

%% ================  Load & prep NODE data 


% % %   loaddir_node           = [locDeriv '/0_finalData_nodes/grp_' parc];
% % % 
% % % % --- Load NODE data
% % %   load([loaddir_node '/CF.mat'],'Dtg');                  CF=Dtg;             % Cortical feature data
% % %   load([loaddir_node '/MR.mat'],'Dtg');                  MR=Dtg;             % Cortical MR data
% % % 
% % % % Labels for cortical data
% % %   lbl_cf={'mean curve' 'medn curve' 'totl curve' 'totl area' 'totl volum' 'mean thick' 'medn thick'};
% % %   lbl_mr={'mean MT' 'medn MT' 'mean R1' 'medn R1'};
% % % 
% % % 
% % % % Info
% % %   Ncf=length(lbl_cf); Nmr=length(lbl_mr);
% % % 
% % % % Trim cortical features
% % %   irm=[1 3 6];
% % %   CF(:,irm)=[]; lbl_cf(irm)=[]; Ncf=length(lbl_cf);
% % % 
% % % % Trim MR
% % %   irm=[1 3];
% % %   MR(:,irm)=[]; lbl_mr(irm)=[]; Nmr=length(lbl_mr);
% % % 
% % % % ----------- Visualize data
% % % % --- Histograms
% % %   ppostmp=getPlotPos([-2559 377 408 420],Ncf,0);
% % % % CF
% % %   D=CF; TMPSTR='Coritcal Features'; TMPLBL=lbl_cf; useN=Ncf;
% % %   for ii=1:useN
% % %     D2=D(:,ii); myfig([TMPSTR strc],ppostmp(ii,:)); histogram(D2,bins); xlabel(TMPLBL{ii}); font(fntsz,usefont);
% % %   end
% % % % MR
% % %   D=MR; TMPSTR='Coritcal MR'; TMPLBL=lbl_mr; useN=Nmr;
% % %   for ii=1:useN
% % %     D2=D(:,ii); myfig([TMPSTR strc],ppostmp(ii,:)); histogram(D2,bins); xlabel(TMPLBL{ii}); font(fntsz,usefont);
% % %   end
% % % % --- Surfaces
% % % % Cortical Features
% % %   D=CF; TMPSTR='Coritcal Features'; TMPLBL=lbl_cf;
% % %   ii=1; D2=D(:,ii); Hcx=plot_conn_surf(D2,pinfo2,'cortex',TMPLBL(ii),TMPSTR); CMAX=round(max(abs(D2)),2); setsurf(Hcx,[CMAX*-1 CMAX],cm2); % Curve
% % %   ii=2; D2=D(:,ii); Hcx=plot_conn_surf(D2,pinfo2,'cortex',TMPLBL(ii),TMPSTR); CMAX=round(max(D2));        setsurf(Hcx,[0 CMAX],cm4);       % Area
% % %   ii=3; D2=D(:,ii); Hcx=plot_conn_surf(D2,pinfo2,'cortex',TMPLBL(ii),TMPSTR); CMAX=round(max(D2));        setsurf(Hcx,[0 CMAX],cm4);       % Volume
% % %   ii=4; D2=D(:,ii); Hcx=plot_conn_surf(D2,pinfo2,'cortex',TMPLBL(ii),TMPSTR); CMAX=round(max(D2),2);      setsurf(Hcx,[0 CMAX],cm4);       % Thick
% % % % Cortical Features
% % %   D=MR; TMPSTR='Coritcal MR'; TMPLBL=lbl_mr;
% % %   ii=1; D2=D(:,ii); Hcx=plot_conn_surf(D2,pinfo2,'cortex',TMPLBL(ii),TMPSTR); CLIM=[min(D2) max(D2)];     setsurf(Hcx,CLIM,cm4);       % MT
% % %   ii=2; D2=D(:,ii); Hcx=plot_conn_surf(D2,pinfo2,'cortex',TMPLBL(ii),TMPSTR); CLIM=[min(D2) max(D2)];     setsurf(Hcx,CLIM,cm4);       % R1
% % % % --- Correlations
% % %   opt={strc,cm3,pposscat,fntsz};
% % % % MR features
% % %   D=MR; TMPLBL=lbl_mr; ddisp_corr(D(:,1),D(:,2),TMPLBL{1},TMPLBL{2},opt{:}); % MT vs R1
% % % % MR features vs S-A Axis
% % %   D=MR; TMPLBL=lbl_mr;
% % %   for ii=1:Nmr; ddisp_corr(SA_map,D(:,ii),'S-A Axis',TMPLBL{ii},opt{:}); end
% % % % Cortical features vs S-A Axis
% % %   D=CF; TMPLBL=lbl_cf;
% % %   for ii=1:Ncf; ddisp_corr(SA_map,D(:,ii),'S-A Axis',TMPLBL{ii},opt{:}); end


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

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%% ============      Communication Models (EDGE data)       ============= %
% -------------------------------------------------------------------------


% Setup
  str1='Communication Models: ';

% ---- Define the data (weighted)
% Other data
  useD={SCnrm MCnrm MC2nrm}; useL={L_sc L_mc L_mc2}; 
  useDstr={'Caliber' extractBefore(MCstr,'-') extractBefore(MC2str,'-')};  useN=length(useL);


% ---- Compute all models for all input data
  ALLCM=cell(1,useN);    CMLBL=cell(1,useN);
  ALLCMnx=cell(1,useN);  CMLBLnx=cell(1,useN);
  for ii=1:useN
      L=useL{ii}; W=useD{ii}; STR=useDstr{ii};
      [ALLCMnx{ii},CMLBLnx{ii},ALLCM{ii},CMLBL{ii}] = communicationModeling(W,L,ED,STR);
  end
  CMLBL_simple=extractAfter(CMLBLnx{1}, '-');
  Ncm=length(CMLBL_simple);

% ---- Repeat for BINARY data
  BIN=double(SC~=0);                                                        % Binary network
  % ALLCM_bin=cell(1,1); CMLBL_bin=cell(1,1); ALLCMnx_bin=cell(1,1); CMLBLnx_bin=cell(1,1);
  [ALLCMnx_bin,CMLBLnx_bin,ALLCM_bin,CMLBL_bin] = communicationModeling_bin(BIN,ED,'Binary');


% ---- Repeat for DELAY data
  % r_delay = 1./delays; r_delay(delays==0)=0;
  [ALLCMnx_delay,CMLBLnx_delay,ALLCM_delay,CMLBL_delay] = communicationModeling(r_delay,delays,ED,delaystr);


%% ------------- BAR PLOTS & SCATTER
% Settings
  % ppostmp=getPlotPos([-2559 489 302 308],Ncm,0); % big
  % ppostmp=getPlotPos([-2559 533 211 264],Ncm,0); % small
  opt4={'XTickLabel',{'all' 'dis' 'con'},'XTickLabelRotation',45};

% % % % ========= Myelin & Caliber communication vs ED for all CMs
% % %   % ppostmp2=getPlotPos([-2559 533 211 264],Ncm,0); % small & thin
% % %   ppostmp2=getPlotPos([-2559 533 250 264],Ncm,0); % wider
% % %   for xx=1:Ncm
% % %       % ---- NO xfm data
% % %       D1=ALLCMnx{1}{xx}; D2=ALLCMnx{2}{xx};
% % %     % -- Correlations with ED
% % %       D3=ED; YSTR='ED'; svcorr=zeros(1,3*2); TMPLBL2=[CMLBL_simple{xx} ' vs ' YSTR];
% % %       X=D1(D3~=0);   Y=D3(D3~=0);                       ii=1; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D2(D3~=0);   Y=D3(D3~=0);                       ii=2; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D1(~masknz & ~maskmd); Y=D3(~masknz & ~maskmd); ii=3; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D2(~masknz & ~maskmd); Y=D3(~masknz & ~maskmd); ii=4; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D1(masknz);  Y=D3(masknz);                      ii=5; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D2(masknz);  Y=D3(masknz);                      ii=6; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %     % Summary bar plot
% % %       myfig([TMPLBL2 str],ppostmp2(xx,:)); B=bar(svcorr,0.6);  hold on; title(CMLBL_simple{xx}); ylim([-.99 .99]);
% % %       ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %       for pp=1:2:6; B.CData(pp,:)=cmap_pred(1,:); end 
% % %       for pp=2:2:6; B.CData(pp,:)=cmap_pred(2,:); end
% % %       B.FaceColor='flat'; set(gca,'XTick',1.5:2:5.5,opt4{:}); 
% % %       for pp=2:2:4; line([pp+.5 pp+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % %       xlim([.2 size(svcorr,2)+.8]);
% % %   end
% % % 
% % % % ========= Myelin, Caliber & DELAY communication vs ED for all CMs
% % %   % ppostmp2=getPlotPos([-2559 533 211 264],Ncm,0); % small & thin
% % %   ppostmp2=getPlotPos([-2559 533 306 264],Ncm,0); % wider
% % %   for xx=1:Ncm
% % %       % ---- NO xfm data
% % %       D1=ALLCMnx{1}{xx}; D2=ALLCMnx{2}{xx}; D4=ALLCMnx_delay{xx};
% % %     % -- Correlations with ED
% % %       D3=ED; YSTR='ED'; svcorr=zeros(1,3*3); TMPLBL2=[CMLBL_simple{xx} ' vs ' YSTR];
% % %       X=D1(D3~=0);   Y=D3(D3~=0);                       ii=1; svcorr(ii)=corr(X,Y,'type','Spearman'); % all
% % %       X=D2(D3~=0);   Y=D3(D3~=0);                       ii=2; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D4(D3~=0);   Y=D3(D3~=0);                       ii=3; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D1(~masknz & ~maskmd); Y=D3(~masknz & ~maskmd); ii=4; svcorr(ii)=corr(X,Y,'type','Spearman'); % dis
% % %       X=D2(~masknz & ~maskmd); Y=D3(~masknz & ~maskmd); ii=5; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D4(~masknz & ~maskmd); Y=D3(~masknz & ~maskmd); ii=6; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D1(masknz);  Y=D3(masknz);                      ii=7; svcorr(ii)=corr(X,Y,'type','Spearman'); % con
% % %       X=D2(masknz);  Y=D3(masknz);                      ii=8; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D4(masknz);  Y=D3(masknz);                      ii=9; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %     % Summary bar plot
% % %       myfig([TMPLBL2 str],ppostmp2(xx,:)); B=bar(svcorr,0.6);  hold on; title(CMLBL_simple{xx}); ylim([-.99 .99]);
% % %       ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %       for pp=1:3:9; B.CData(pp,:)=cmap_pred(1,:); end 
% % %       for pp=2:3:9; B.CData(pp,:)=cmap_pred(2,:); end
% % %       for pp=3:3:9; B.CData(pp,:)=cmap_pred(4,:); end
% % %       B.FaceColor='flat'; set(gca,'XTick',2:3:8,opt4{:}); 
% % %       for pp=3:3:6; line([pp+.5 pp+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % %       xlim([.2 size(svcorr,2)+.8]);
% % %   end

% ======== Myelin communication vs Caliber
% Bar plot correlations between communication networks
  ppostmp=getPlotPos([-2559 577 164 220],Ncm,0); % smallest
  for xx=1:Ncm
      % ---- NO xfm data
      D1=ALLCMnx{1}{xx};  Dstr1=CMLBLnx{1}{xx};
      D2=ALLCMnx{2}{xx};  Dstr2=CMLBLnx{2}{xx};
      % Db=ALLCMnx_bin{xx}; Dstrb=CMLBLnx_bin{xx}; CLIMb=[prctile(Db(Db~=0),1)  prctile(Db(Db~=0),99)];
    % -- Correlations with caliber
      svcorr=zeros(1,3); TMPLBL2=['caliber vs myelin (' CMLBL_simple{xx} ')'];
      X=D1(~maskmd);           Y=D2(~maskmd);           ii=1; svcorr(ii)=corr(X,Y,'type','Spearman');
      X=D1(~masknz & ~maskmd); Y=D2(~masknz & ~maskmd); ii=2; svcorr(ii)=corr(X,Y,'type','Spearman');
      X=D1(masknz);            Y=D2(masknz);            ii=3; svcorr(ii)=corr(X,Y,'type','Spearman');
    % Summary bar plot
      myfig([TMPLBL2 str],ppostmp(xx,:)); B=bar(svcorr,0.6);  hold on; ylim([-.99 .99])
      font(fntsz-2,usefont); grid on; set(gca,'XGrid','off'); % ylabel('\rho');
      B.CData=cmap_6cm(xx,:); B.FaceColor='flat'; set(gca,opt4{:}); % title(CMLBL_simple{xx});
      xlim([.3 size(svcorr,2)+.7]); set(gca,'TickLength',[0 0]);
  end


% ======== DELAY communication vs Caliber
  ppostmp=getPlotPos([-2559 577 164 220],Ncm,0); % smallest
  for xx=1:Ncm
      % ---- NO xfm data
      D1=ALLCMnx{1}{xx};     Dstr1=CMLBLnx{1}{xx};
      D2=ALLCMnx_delay{xx};  Dstr2=CMLBLnx_delay{xx};
      % Db=ALLCMnx_bin{xx}; Dstrb=CMLBLnx_bin{xx}; CLIMb=[prctile(Db(Db~=0),1)  prctile(Db(Db~=0),99)];
    % -- Correlations with caliber
      svcorr=zeros(1,3); TMPLBL2=['caliber vs delay (' CMLBL_simple{xx} ')'];
      X=D1(~maskmd);           Y=D2(~maskmd);           ii=1; svcorr(ii)=corr(X,Y,'type','Spearman');
      X=D1(~masknz & ~maskmd); Y=D2(~masknz & ~maskmd); ii=2; svcorr(ii)=corr(X,Y,'type','Spearman');
      X=D1(masknz);            Y=D2(masknz);            ii=3; svcorr(ii)=corr(X,Y,'type','Spearman');
    % Summary bar plot
      myfig([TMPLBL2 str],ppostmp(xx,:)); B=bar(svcorr,0.6);  hold on; ylim([-.99 .99])
      font(fntsz-2,usefont); grid on; set(gca,'XGrid','off'); % ylabel('\rho');
      B.CData=cmap_6cm(xx,:); B.FaceColor='flat'; set(gca,opt4{:}); % title(CMLBL_simple{xx});
      xlim([.3 size(svcorr,2)+.7]); set(gca,'TickLength',[0 0]); nolbl;
  end

% % % % ======== Myelin communication vs Binary
% % %   ppostmp=getPlotPos([-2559 577 164 220],Ncm,0); % smallest
% % %   for xx=1:Ncm
% % %       % ---- NO xfm data
% % %       D2=ALLCMnx{2}{xx}; D1=ALLCMnx_bin{xx}; 
% % %     % -- Correlations with caliber
% % %       svcorr=zeros(1,3); TMPLBL2='BINARY vs myelin';
% % %       X=D1(~maskmd);           Y=D2(~maskmd);           ii=1; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D1(~masknz & ~maskmd); Y=D2(~masknz & ~maskmd); ii=2; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %       X=D1(masknz);            Y=D2(masknz);            ii=3; svcorr(ii)=corr(X,Y,'type','Spearman');
% % %     % Summary bar plot
% % %       myfig([TMPLBL2 str],ppostmp(xx,:)); B=bar(svcorr,0.6);  hold on; title(CMLBL_simple{xx}); ylim([-1 1]);
% % %       ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %       B.CData=cmap_6cm(xx,:); B.FaceColor='flat'; set(gca,opt4{:}); 
% % %       xlim([.3 size(svcorr,2)+.7]);
% % %   end
% % % 
% % % % ==========   Scatter plot caliber vs MYELIN: data colored by connectivity
% % %   tmask=zeros(Nnode); tmask(masknz)=1; ppostmp=getPlotPos([-2559 489 302 308],Ncm,0);
% % %   TMPLBL2='CALIBER vs myelin';
% % %   for xx=1:Ncm 
% % %     % ---- zscored data
% % %       D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx};
% % %     % -- All edges colored by connectivity
% % %       X=D1(:); Y=D2(:); tmask_vec=tmask(:);                                     % flatten
% % %       X(maskmd(:))=[]; Y(maskmd(:))=[]; tmask_vec(maskmd(:))=[];                % Remove main diagonal
% % %     % Function to compute density and scatter
% % %       plotDensityScatt=@(x,y,baseColor) scatter(x,y,5,densityColor(x,y,baseColor),'filled','MarkerFaceAlpha',0.8);
% % %     % Plot
% % %       myfig(TMPLBL2,ppostmp(xx,:)); hold on;
% % %       plotDensityScatt(X(tmask_vec==1), Y(tmask_vec==1), cmap_discon(2,:)); % red = connected
% % %       plotDensityScatt(X(~tmask_vec),Y(~tmask_vec),cmap_discon(1,:)); % blue = disconnected
% % %       xlabel(extractBefore(Dstr1,'-')); ylabel(extractBefore(Dstr2,'-')); title(CMLBL_simple{xx}); font(fntsz,usefont);
% % %       box on;
% % %   end




% % % %% ------------- SCATTER PLOTS + HISTOGRAMS
% % % % Settings
% % %   % ppostmp=getPlotPos([-2559 489 302 308],Ncm,0); % large
% % %   ppostmp=getPlotPos([-2559 562 243 235],Ncm,0); % small
% % %   OPTX={'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability'};
% % %   OPTY={'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability','Orientation','horizontal'};
% % %   Ntmp=4;  % governs ratio of histograms to scatter
% % % 
% % % % ==========   Scatter plot CALIBER vs myelin: data colored by connectivity (+ histograms)
% % %   tmask=zeros(Nnode); tmask(masknz & maskut)=1; tmask(~maskut)=-1; tmask(maskmd)=-1; 
% % %   XLIMS={[-2.5 5], [-5 5],     [-2.2 6.2], [-2.7 5],   [-3.2 2.9], [-3.6 10.6]};
% % %   YLIMS={[-2.5 5], [-4.9 4.9], [-1.9 3.9], [-3.9 5.9], [-4 2.5],   [-3.9 4.2]};
% % %   for xx=1:Ncm 
% % %     % ---- zscored data
% % %       D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx};
% % %     % -- All edges colored by connectivity
% % %       Xc=D1(tmask==1); Xd=D1(tmask==0); Yc=D2(tmask==1); Yd=D2(tmask==0);    % flatten
% % %     % Function to compute density and scatter
% % %       plotDensityScatt=@(x,y,baseColor) scatter(x,y,5,densityColor(x,y,baseColor),'filled','MarkerFaceAlpha',0.8);
% % %     % --- Plot
% % %       myfig(['CALIBER vs myelin: ' CMLBL_simple{xx}],ppostmp(xx,:));
% % %       tl=tiledlayout(Ntmp,Ntmp,'TileSpacing','compact','Padding','compact'); % Tiled layout for scatter + marginals
% % %     % Histograms on top (X-axis)
% % %       axHistX=nexttile(tl,1,[1 Ntmp-1]); hold on;
% % %       histogram(axHistX,Xc,'FaceColor',cmap_discon(2,:),OPTX{:});
% % %       histogram(axHistX,Xd,'FaceColor',cmap_discon(1,:),OPTX{:});
% % %       axHistX.XAxisLocation='top'; axHistX.XAxis.Visible='off'; axHistX.YAxis.Visible='off';
% % %     % Scatter (main panel)
% % %       axScatter=nexttile(tl,Ntmp+1,[Ntmp-1 Ntmp-1]); hold on;
% % %       plotDensityScatt(Xc, Yc, cmap_discon(2,:));                           % red = connected
% % %       plotDensityScatt(Xd,Yd,cmap_discon(1,:));                             % blue = disconnected
% % %       % xlabel(axScatter,extractBefore(Dstr1,'-')); ylabel(axScatter,extractBefore(Dstr2,'-')); font(fntsz-2,usefont); %title(CMLBL_simple{xx});
% % %       box(axScatter,'on'); xlim(XLIMS{xx}); ylim(YLIMS{xx}); 
% % %       nolbl
% % %     % Histograms on right (Y-axis)
% % %       axHistY=nexttile(tl,Ntmp*2,[Ntmp-1 1]); hold on;
% % %       histogram(axHistY, Yc,'FaceColor',cmap_discon(2,:),OPTY{:});
% % %       histogram(axHistY, Yd,'FaceColor',cmap_discon(1,:),OPTY{:});
% % %       axHistY.YAxisLocation='right'; axHistY.XAxis.Visible='off'; axHistY.YAxis.Visible='off';
% % %     % Link axes & line up lims
% % %       axHistX.XLim=axScatter.XLim; axHistY.YLim=axScatter.YLim;
% % %       linkaxes([axScatter,axHistX],'x'); linkaxes([axScatter,axHistY],'y');
% % %   end
% % % % TMP SAVE
% % %   set(gcf,'Renderer','painters'); svstr=[num2str(pp) '_' CMLBL_simple{pp}];
% % %   exportgraphics(gcf,['~/Downloads/3-mtsat-scat_' svstr '_nolbl.pdf'],'ContentType','vector');
% % % 
% % % 
% % % % ==========   Scatter plot CALIBER vs DELAY: data colored by connectivity (+ histograms)
% % %   tmask=zeros(Nnode); tmask(masknz & maskut)=1; tmask(~maskut)=-1; tmask(maskmd)=-1; 
% % %   XLIMS={[-2.1 7], [-4 7],   [-2.2 5.9], [-2.7 5], [-3.2 2.9], [-3.6 7.9]};
% % %   YLIMS={[-2 5.5], [-3.5 5], [-2.3 3.5], [-3 5],   [-4 3],     [-3.9 5.9]};
% % %   for xx=1:Ncm 
% % %     % ---- zscored data
% % %       D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; D2=ALLCM_delay{xx}; Dstr2=CMLBL_delay{xx};
% % %     % -- All edges colored by connectivity
% % %       Xc=D1(tmask==1); Xd=D1(tmask==0); Yc=D2(tmask==1); Yd=D2(tmask==0);    % flatten
% % %     % Function to compute density and scatter
% % %       plotDensityScatt=@(x,y,baseColor) scatter(x,y,5,densityColor(x,y,baseColor),'filled','MarkerFaceAlpha',0.8);
% % %     % --- Plot
% % %       myfig(['CALIBER vs myelin: ' CMLBL_simple{xx}],ppostmp(xx,:));
% % %     % Tiled layout for scatter + marginals
% % %       tl = tiledlayout(Ntmp,Ntmp,'TileSpacing','compact','Padding','compact');
% % %     % Histograms on top (X-axis)
% % %       axHistX=nexttile(tl,1,[1 Ntmp-1]); hold on;
% % %       histogram(axHistX,Xc,'FaceColor',cmap_discon(2,:),OPTX{:});
% % %       histogram(axHistX,Xd,'FaceColor',cmap_discon(1,:),OPTX{:});
% % %       axHistX.XAxisLocation='top'; axHistX.XAxis.Visible='off'; axHistX.YAxis.Visible='off';
% % %     % Scatter (main panel)
% % %       axScatter=nexttile(tl,Ntmp+1,[Ntmp-1 Ntmp-1]); hold on;
% % %       plotDensityScatt(Xc, Yc, cmap_discon(2,:));                           % red = connected
% % %       plotDensityScatt(Xd,Yd,cmap_discon(1,:));                             % blue = disconnected
% % %       % xlabel(axScatter,extractBefore(Dstr1,'-')); ylabel(axScatter,extractBefore(Dstr2,'-')); font(fntsz-2,usefont); %title(CMLBL_simple{xx});
% % %       box(axScatter,'on'); xlim(XLIMS{xx}); ylim(YLIMS{xx});
% % %       nolbl
% % %     % Histograms on right (Y-axis)
% % %       axHistY=nexttile(tl,Ntmp*2,[Ntmp-1 1]); hold on;
% % %       histogram(axHistY, Yc,'FaceColor',cmap_discon(2,:),OPTY{:});
% % %       histogram(axHistY, Yd,'FaceColor',cmap_discon(1,:),OPTY{:});
% % %       axHistY.YAxisLocation='right'; axHistY.XAxis.Visible='off'; axHistY.YAxis.Visible='off';
% % %     % Link axes & line up lims
% % %       axHistX.XLim=axScatter.XLim; axHistY.YLim=axScatter.YLim;
% % %       linkaxes([axScatter,axHistX],'x'); linkaxes([axScatter,axHistY],'y');
% % %   end



  %% ------------- SCATTER PLOTS + HISTOGRAMS (BINNED DENSITY)
% Settings
  % ppostmp=getPlotPos([-2559 489 302 308],Ncm,0); % large
  ppostmp=getPlotPos([-2559 562 243 235],Ncm,0); % small
  OPTX={'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability'};
  OPTY={'FaceAlpha',0.5,'EdgeColor','none','Normalization','probability','Orientation','horizontal'};
  Ntmp=4;  % governs ratio of histograms to scatter

% ==========   Scatter plot CALIBER vs myelin: data colored by connectivity (+ histograms)
  tmask=zeros(Nnode); tmask(masknz & maskut)=1; tmask(~maskut)=-1; tmask(maskmd)=-1; 
  XLIMS={[-2.5 5], [-5 5],     [-2.2 6.2], [-2.7 5],   [-3.2 2.9], [-3.6 10.6]};
  YLIMS={[-2.5 5], [-4.9 4.9], [-1.9 3.9], [-3.9 5.9], [-4 2.5],   [-3.9 4.2]};
  densopt = struct;
  densopt.nbins = 200;         % increase for smoother / more detailed
  densopt.sigma = 1.0;         % gaussian smoothing in bin units (0 = none)
  densopt.alphaMax = 0.95;     % cap opacity
  densopt.gamma = 0.9;         % <1 boosts low density visibility
  for xx=1:Ncm 
    % ---- zscored data
      D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx};
    % -- All edges colored by connectivity
      Xc=D1(tmask==1); Xd=D1(tmask==0); Yc=D2(tmask==1); Yd=D2(tmask==0);    % flatten
    % --- Plot
      myfig(['CALIBER vs myelin: ' CMLBL_simple{xx}],ppostmp(xx,:));
      tl=tiledlayout(Ntmp,Ntmp,'TileSpacing','compact','Padding','compact'); % Tiled layout for scatter + marginals
    % Histograms on top (X-axis)
      axHistX=nexttile(tl,1,[1 Ntmp-1]); hold on;
      histogram(axHistX,Xc,'FaceColor',cmap_discon(2,:),OPTX{:});
      histogram(axHistX,Xd,'FaceColor',cmap_discon(1,:),OPTX{:});
      axHistX.XAxisLocation='top'; axHistX.XAxis.Visible='off'; axHistX.YAxis.Visible='off';
    % Scatter (main panel)
      axScatter=nexttile(tl,Ntmp+1,[Ntmp-1 Ntmp-1]); hold on;
      densopt.xlim = XLIMS{xx};  densopt.ylim = YLIMS{xx};
      plotDensity2D_overlay(axScatter, Xd, Yd, cmap_discon(1,:), densopt); % blue = disconnected
      plotDensity2D_overlay(axScatter, Xc, Yc, cmap_discon(2,:), densopt); % red = connected
      % xlabel(axScatter,extractBefore(Dstr1,'-')); ylabel(axScatter,extractBefore(Dstr2,'-')); font(fntsz-2,usefont); %title(CMLBL_simple{xx});
      box(axScatter,'on'); xlim(XLIMS{xx}); ylim(YLIMS{xx}); 
      nolbl
    % Histograms on right (Y-axis)
      axHistY=nexttile(tl,Ntmp*2,[Ntmp-1 1]); hold on;
      histogram(axHistY, Yc,'FaceColor',cmap_discon(2,:),OPTY{:});
      histogram(axHistY, Yd,'FaceColor',cmap_discon(1,:),OPTY{:});
      axHistY.YAxisLocation='right'; axHistY.XAxis.Visible='off'; axHistY.YAxis.Visible='off';
    % Link axes & line up lims
      axHistX.XLim=axScatter.XLim; axHistY.YLim=axScatter.YLim;
      linkaxes([axScatter,axHistX],'x'); linkaxes([axScatter,axHistY],'y');
  end
% TMP SAVE
  set(gcf,'Renderer','painters'); svstr=[num2str(pp) '_' CMLBL_simple{pp}];
  exportgraphics(gcf,['~/Downloads/3-mtsat-scat_' svstr '_nolbl.pdf'],'ContentType','vector');


% ==========   Scatter plot CALIBER vs DELAY: data colored by connectivity (+ histograms)
  tmask=zeros(Nnode); tmask(masknz & maskut)=1; tmask(~maskut)=-1; tmask(maskmd)=-1; 
  XLIMS={[-2.1 7], [-4 7],   [-2.2 5.9], [-2.7 5], [-3.2 2.9], [-3.6 7.9]};
  YLIMS={[-2 5.5], [-3.5 5], [-2.3 3.5], [-3 5],   [-4 3],     [-3.9 5.9]};
  for xx=1:Ncm 
    % ---- zscored data
      D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; D2=ALLCM_delay{xx}; Dstr2=CMLBL_delay{xx};
    % -- All edges colored by connectivity
      Xc=D1(tmask==1); Xd=D1(tmask==0); Yc=D2(tmask==1); Yd=D2(tmask==0);    % flatten
    % --- Plot
      myfig(['CALIBER vs myelin: ' CMLBL_simple{xx}],ppostmp(xx,:));
    % Tiled layout for scatter + marginals
      tl = tiledlayout(Ntmp,Ntmp,'TileSpacing','compact','Padding','compact');
    % Histograms on top (X-axis)
      axHistX=nexttile(tl,1,[1 Ntmp-1]); hold on;
      histogram(axHistX,Xc,'FaceColor',cmap_discon(2,:),OPTX{:});
      histogram(axHistX,Xd,'FaceColor',cmap_discon(1,:),OPTX{:});
      axHistX.XAxisLocation='top'; axHistX.XAxis.Visible='off'; axHistX.YAxis.Visible='off';
    % Scatter (main panel)
      axScatter=nexttile(tl,Ntmp+1,[Ntmp-1 Ntmp-1]); hold on;
      densopt.xlim = XLIMS{xx};  densopt.ylim = YLIMS{xx};
      plotDensity2D_overlay(axScatter, Xd, Yd, cmap_discon(1,:), densopt); % blue = disconnected
      plotDensity2D_overlay(axScatter, Xc, Yc, cmap_discon(2,:), densopt); % red = connected
      % xlabel(axScatter,extractBefore(Dstr1,'-')); ylabel(axScatter,extractBefore(Dstr2,'-')); font(fntsz-2,usefont); %title(CMLBL_simple{xx});
      box(axScatter,'on'); xlim(XLIMS{xx}); ylim(YLIMS{xx});
      nolbl
    % Histograms on right (Y-axis)
      axHistY=nexttile(tl,Ntmp*2,[Ntmp-1 1]); hold on;
      histogram(axHistY, Yc,'FaceColor',cmap_discon(2,:),OPTY{:});
      histogram(axHistY, Yd,'FaceColor',cmap_discon(1,:),OPTY{:});
      axHistY.YAxisLocation='right'; axHistY.XAxis.Visible='off'; axHistY.YAxis.Visible='off';
    % Link axes & line up lims
      axHistX.XLim=axScatter.XLim; axHistY.YLim=axScatter.YLim;
      linkaxes([axScatter,axHistX],'x'); linkaxes([axScatter,axHistY],'y');
  end
% TMP SAVE
  set(gcf,'Renderer','painters'); svstr=[num2str(pp) '_' CMLBL_simple{pp}];
  exportgraphics(gcf,['~/Downloads/4_delay-scat_' svstr '_nolbl.pdf'],'ContentType','vector');



% % % % ==========   Scatter plot CALIBER vs myelin: data colored by ED (+ histograms)
% % %   XLIMS={[-2.5 5], [-5 5], [-2.2 6.2], [-2.7 5], [-3.2 2.9], [-3.6 10.6]};
% % %   YLIMS={[-2.5 5], [-5 5], [-2 4],     [-4 6],   [-4 2.3],   [-4 4.2]};
% % %   Ntmp=4;  D3=nzzscore(EDnrm,0); 
% % %   mask1=D3>0; mask2=D3<0; mask1=mask1(~maskmd(:)); mask2=mask2(~maskmd(:));
% % %   for xx=1:Ncm 
% % %     % ---- zscored data
% % %       D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx};
% % %     % -- All edges colored by ED
% % %       X=D1(~maskmd(:)); Y=D2(~maskmd(:)); Z=D3(~maskmd(:));                                                % flatten
% % %       % X(maskmd(:))=[]; Y(maskmd(:))=[]; Z(maskmd(:))=[];                        % Remove main diagonal
% % %     % Function to compute density and scatter
% % %       plotDensityScatt=@(x,y,baseColor) scatter(x,y,5,densityColor(x,y,baseColor),'filled','MarkerFaceAlpha',0.8);
% % %     % --- Plot
% % %       myfig(['CALIBER vs myelin (ED): ' CMLBL_simple{xx}],ppostmp(xx,:));
% % %     % Tiled layout for scatter + marginals
% % %       tl = tiledlayout(Ntmp,Ntmp,'TileSpacing','compact','Padding','compact');
% % %     % Histograms on top (X-axis)
% % %       axHistX=nexttile(tl,1,[1 Ntmp-1]); hold on;
% % %       histogram(axHistX,X(mask2),'FaceColor',cmap_discon(1,:),OPTX{:});      % blue = ED < mean
% % %       histogram(axHistX,X(mask1),'FaceColor',cmap_discon(2,:),OPTX{:});      % red  = ED > mean
% % %       axHistX.XAxisLocation='top'; axHistX.XAxis.Visible='off'; axHistX.YAxis.Visible='off';
% % %     % Scatter (main panel)
% % %       axScatter=nexttile(tl,Ntmp+1,[Ntmp-1 Ntmp-1]); hold on;
% % %       plotDensityScatt(X(mask2),Y(mask2),cmap_discon(1,:));                  % blue = ED < mean
% % %       plotDensityScatt(X(mask1),Y(mask1),cmap_discon(2,:));                  % red  = ED > mean
% % %       xlabel(axScatter,extractBefore(Dstr1,'-')); ylabel(axScatter,extractBefore(Dstr2,'-')); font(fntsz,usefont); %title(CMLBL_simple{xx});
% % %       box(axScatter,'on'); xlim(XLIMS{xx}); ylim(YLIMS{xx});
% % %       line(xlim,[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth',1.5);
% % %       line([0 0],ylim,'Color',[0 0 0],'LineStyle','--','LineWidth',1.5);
% % %     % Histograms on right (Y-axis)
% % %       axHistY=nexttile(tl,Ntmp*2,[Ntmp-1 1]); hold on;
% % %       histogram(axHistY,Y(mask2),'FaceColor',cmap_discon(1,:),OPTY{:});      % blue = ED < mean
% % %       histogram(axHistY,Y(mask1),'FaceColor',cmap_discon(2,:),OPTY{:});      % red  = ED > mean
% % %       axHistY.YAxisLocation='right'; axHistY.XAxis.Visible='off'; axHistY.YAxis.Visible='off';
% % %     % Link axes & line up lims
% % %       axHistX.XLim=axScatter.XLim; axHistY.YLim=axScatter.YLim;
% % %       linkaxes([axScatter,axHistX],'x'); linkaxes([axScatter,axHistY],'y');
% % %   end


%% ------------- BOX PLOTS: comparison to distance
% % % % ==========   BOX plot CALIBER vs myelin: ED separated by < & > MEAN 
% % %   ppostmp=getPlotPos([-2559 443 305 354],Ncm,0); ppostmp2=getPlotPos([-2559 9 305 354],Ncm,0);
% % %   scatopt={'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','w','LineWidth',0.1};
% % %   TLTOPT={'FontWeight','normal','FontSize',20};
% % %   D3=nzzscore(EDnrm,0); mask1=D3>0; mask2=D3<0; mask1=mask1(~maskmd(:)); mask2=mask2(~maskmd(:));
% % %   for xx=1:Ncm 
% % %       Dt1=ALLCM{1}{xx}; Dt2=ALLCM{2}{xx};                                   % zscored data
% % %       X=Dt1(~maskmd(:)); Y=Dt2(~maskmd(:));                                 % All edges excluding main diagonal
% % % 
% % %     % ----- ED > mean
% % %       % tstats=zeros(2,Ncis); cohens_d=zeros(2,Ncis);
% % %     % Prep
% % %       D1=X(mask1);   D2=Y(mask1);
% % %       idx=unique([find(isnan(D1))   find(isnan(D2))]);
% % %       D1(idx)=[];          D2(idx)=[];
% % %       DPLT=[D1; D2];                                                        % Combine data into one vector
% % %       tlbls=[]; for xxx=1:2; tlbls=[tlbls; xxx*ones(length(D1),1)]; end       % Create group labels
% % %     % t-test
% % %       [h,p,ci,stats]=ttest2(D1,D2);
% % %     % Cohen's d
% % %       t_m1=mean(D1);   t_m2=mean(D2);
% % %       t_s1=std(D1);    t_s2=std(D2);
% % %       t_n1=length(D1); t_n2=length(D2);
% % %       t_pooled_sd=sqrt(((t_n1-1)*t_s1^2+(t_n2-1)*t_s2^2)/(t_n1+t_n2-2));
% % %       t_cohens_d=(t_m1-t_m2)/t_pooled_sd;
% % %     % Plot
% % %       myfig(['ED > mean: ' CMLBL_simple{xx}],ppostmp(xx,:)); hold on; colormap(cm2); clim([-3 3]);
% % %       x_jit=0.05; xpos=cell(1,2); for xxx=1:2; xpos{xxx}=xxx+x_jit*randn(length(D1),1); end  % Jitter x positions for better visibility
% % %     % Scatter plot with color-coded residuals
% % %       scatter(xpos{1},D1,10,D1,scatopt{:});
% % %       scatter(xpos{2},D2,10,D2,scatopt{:}); % colorbar;
% % %     % Add boxplots on top (no symbols)
% % %       B=boxplot(DPLT,tlbls,'Colors','k','Widths',0.4,'Symbol','');
% % %       set(findobj(B,'Type','Line'),'LineWidth', 2); font(fntsz,usefont);
% % %       set(gca,'XTick',1:2,'XTickLabel',{'Caliber' 'Myelin'}); ylabel(CMLBL_simple{xx}); box on;
% % %       grid on; set(gca,'XGrid','off');
% % %       title([sprintf('t=%.1f, d=%.2f',stats.tstat,t_cohens_d)],TLTOPT{:});
% % %     % ----- ED < mean
% % %       % tstats=zeros(2,Ncis); cohens_d=zeros(2,Ncis);
% % %     % Prep
% % %       D1=X(mask2);   D2=Y(mask2);
% % %       idx=unique([find(isnan(D1))   find(isnan(D2))]);
% % %       D1(idx)=[];          D2(idx)=[];
% % %       DPLT=[D1; D2];                                                        % Combine data into one vector
% % %       tlbls=[]; for xxx=1:2; tlbls=[tlbls; xxx*ones(length(D1),1)]; end       % Create group labels
% % %     % t-test
% % %       [h,p,ci,stats]=ttest2(D1,D2);
% % %     % Cohen's d
% % %       t_m1=mean(D1);   t_m2=mean(D2);
% % %       t_s1=std(D1);    t_s2=std(D2);
% % %       t_n1=length(D1); t_n2=length(D2);
% % %       t_pooled_sd=sqrt(((t_n1-1)*t_s1^2+(t_n2-1)*t_s2^2)/(t_n1+t_n2-2));
% % %       t_cohens_d=(t_m1-t_m2)/t_pooled_sd;
% % %     % Plot
% % %       myfig(['ED < mean: ' CMLBL_simple{xx}],ppostmp2(xx,:)); hold on; colormap(cm2); clim([-3 3]);
% % %       x_jit=0.05; xpos=cell(1,2); for xxx=1:2; xpos{xxx}=xxx+x_jit*randn(length(D1),1); end  % Jitter x positions for better visibility
% % %     % Scatter plot with color-coded residuals
% % %       scatter(xpos{1},D1,10,D1,scatopt{:});
% % %       scatter(xpos{2},D2,10,D2,scatopt{:}); % colorbar;
% % %     % Add boxplots on top (no symbols)
% % %       B=boxplot(DPLT,tlbls,'Colors','k','Widths',0.4,'Symbol','');
% % %       set(findobj(B,'Type','Line'),'LineWidth', 2); font(fntsz,usefont);
% % %       set(gca,'XTick',1:2,'XTickLabel',{'Caliber' 'Myelin'}); ylabel(CMLBL_simple{xx}); box on;
% % %       grid on; set(gca,'XGrid','off');
% % %       title([sprintf('t=%.1f, d=%.2f',stats.tstat,t_cohens_d)],TLTOPT{:});
% % %   end

% % % % ==========   BOX plot CALIBER vs myelin: ED separated by 1st & 3rd quartiles 
% % % % Switched order of caliber & myelin
% % %   % ppostmp=getPlotPos([-2559 443 305 354],Ncm,0); ppostmp2=getPlotPos([-2559 9 305 354],Ncm,0); % big
% % %   % ppostmp=getPlotPos([-2559 496 272 301],Ncm,0); ppostmp2=getPlotPos([-2559 115 272 301],Ncm,0); % small
% % %   ppostmp=getPlotPos([-2559 548 230 249],Ncm,0); ppostmp2=getPlotPos([-2559 185 230 249],Ncm,0); % smaller
% % %   scatopt={'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','w','LineWidth',0.1};
% % %   TLTOPT={'FontWeight','normal','FontSize',20};
% % %   D3=EDnrm; Q1=prctile(D3(D3~=0),25); Q3=prctile(D3(D3~=0),75);
% % %   mask1=D3>Q3; mask2=D3<Q1; mask1=mask1(~maskmd(:)); mask2=mask2(~maskmd(:));
% % %   tstats=zeros(2,Ncm); cohens_d=zeros(2,Ncm);
% % %   for xx=1:Ncm 
% % %       Dt1=ALLCM{2}{xx}; Dt2=ALLCM{1}{xx};                                   % zscored data
% % %       X=Dt1(~maskmd(:)); Y=Dt2(~maskmd(:));                                 % All edges excluding main diagonal
% % % 
% % %     % ----- ED > Q3
% % %     % Prep
% % %       D1=X(mask1);   D2=Y(mask1);
% % %       idx=unique([find(isnan(D1))   find(isnan(D2))]);
% % %       D1(idx)=[];          D2(idx)=[];
% % %       DPLT=[D1; D2];                                                        % Combine data into one vector
% % %       tlbls=[]; for xxx=1:2; tlbls=[tlbls; xxx*ones(length(D1),1)]; end       % Create group labels
% % %     % t-test
% % %       [h,p,ci,stats]=ttest2(D1,D2);
% % %     % Cohen's d
% % %       t_m1=mean(D1);   t_m2=mean(D2);
% % %       t_s1=std(D1);    t_s2=std(D2);
% % %       t_n1=length(D1); t_n2=length(D2);
% % %       t_pooled_sd=sqrt(((t_n1-1)*t_s1^2+(t_n2-1)*t_s2^2)/(t_n1+t_n2-2));
% % %       t_cohens_d=(t_m1-t_m2)/t_pooled_sd;
% % %     % Plot
% % %       myfig(['ED > Q3: ' CMLBL_simple{xx}],ppostmp(xx,:)); hold on; colormap(cm2); clim([-3 3]);
% % %       yline([0 0],'Color',[0.5 0.5 0.5],'LineWidth',1);
% % %       x_jit=0.05; xpos=cell(1,2); for xxx=1:2; xpos{xxx}=xxx+x_jit*randn(length(D1),1); end  % Jitter x positions for better visibility
% % %     % Scatter plot with color-coded residuals
% % %       scatter(xpos{1},D1,10,D1,scatopt{:});
% % %       scatter(xpos{2},D2,10,D2,scatopt{:}); % colorbar;
% % %     % Add boxplots on top (no symbols)
% % %       B=boxplot(DPLT,tlbls,'Colors','k','Widths',0.4,'Symbol','');
% % %       set(findobj(B,'Type','Line'),'LineWidth', 2); font(fntsz,usefont);
% % %       set(gca,'XTick',1:2,'XTickLabel',{'Myelin' 'Caliber'},'XTickLabelRotation',0); 
% % %       ylabel(CMLBL_simple{xx}); box on;
% % %       grid on; set(gca,'XGrid','off'); set(gca,'TickLength',[0 0]);
% % %       % title([sprintf('t=%.1f, d=%.2f',stats.tstat,t_cohens_d)],TLTOPT{:});
% % %       title([sprintf('d = %.2f',t_cohens_d)],TLTOPT{:});      
% % %     % Save stats
% % %       tstats(1,xx)=stats.tstat; cohens_d(1,xx)=t_cohens_d;
% % %     % ----- ED < Q1
% % %     % Prep
% % %       D1=X(mask2);   D2=Y(mask2);
% % %       idx=unique([find(isnan(D1))   find(isnan(D2))]);
% % %       D1(idx)=[];          D2(idx)=[];
% % %       DPLT=[D1; D2];                                                        % Combine data into one vector
% % %       tlbls=[]; for xxx=1:2; tlbls=[tlbls; xxx*ones(length(D1),1)]; end       % Create group labels
% % %     % t-test
% % %       [h,p,ci,stats]=ttest2(D1,D2);
% % %     % Cohen's d
% % %       t_m1=mean(D1);   t_m2=mean(D2);
% % %       t_s1=std(D1);    t_s2=std(D2);
% % %       t_n1=length(D1); t_n2=length(D2);
% % %       t_pooled_sd=sqrt(((t_n1-1)*t_s1^2+(t_n2-1)*t_s2^2)/(t_n1+t_n2-2));
% % %       t_cohens_d=(t_m1-t_m2)/t_pooled_sd;
% % %     % Plot
% % %       myfig(['ED < Q1: ' CMLBL_simple{xx}],ppostmp2(xx,:)); hold on; colormap(cm2); clim([-3 3]);
% % %       yline([0 0],'Color',[0.5 0.5 0.5],'LineWidth',1);
% % %       x_jit=0.05; xpos=cell(1,2); for xxx=1:2; xpos{xxx}=xxx+x_jit*randn(length(D1),1); end  % Jitter x positions for better visibility
% % %     % Scatter plot with color-coded residuals
% % %       scatter(xpos{1},D1,10,D1,scatopt{:});
% % %       scatter(xpos{2},D2,10,D2,scatopt{:}); % colorbar;
% % %     % Add boxplots on top (no symbols)
% % %       B=boxplot(DPLT,tlbls,'Colors','k','Widths',0.4,'Symbol','');
% % %       set(findobj(B,'Type','Line'),'LineWidth', 2); font(fntsz,usefont);
% % %       set(gca,'XTick',1:2,'XTickLabel',{'Myelin' 'Caliber'},'XTickLabelRotation',0); 
% % %       ylabel(CMLBL_simple{xx}); box on;
% % %       grid on; set(gca,'XGrid','off'); set(gca,'TickLength',[0 0]);
% % %       % title([sprintf('t=%.1f, d=%.2f',stats.tstat,t_cohens_d)],TLTOPT{:});
% % %       title([sprintf('d = %.2f',t_cohens_d)],TLTOPT{:});      
% % %     % Save stats
% % %       tstats(2,xx)=stats.tstat; cohens_d(2,xx)=t_cohens_d;
% % %   end
% % % % Summary: t-stats
% % %   D=tstats; ppostmp=[-2559 600 569 197]; % ppostmp=[-2559 537 687 260];
% % %   YLIM=[floor(min(D(:)))  ceil(max(D(:)))]; % YLIM=[floor(min(D(:))*10)/10  ceil(max(D(:))*10)/10];
% % %   myfig('Communication Efficiency vs dist',ppostmp); B=bar(D',0.8); hold on; font(fntsz,usefont);
% % %   ylabel('t-stat'); grid on; set(gca,'XGrid','off'); ylim(YLIM); 
% % %   B(1).FaceColor=[0 0 0]; B(1).FaceAlpha=1; B(2).FaceColor=[0 0 0]; B(2).FaceAlpha=.3;
% % %   for ii=1:Ncm-1; line([ii+.5 ii+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % %   set(gca,'XTickLabels',CMLBL_simple); 
% % % % Summary: Cohens d
% % %   D=cohens_d; ppostmp=[-2559 600 569 197]; %ppostmp=[-2559 537 687 260];
% % %   YLIM=[-1 1]; %YLIM=[floor(min(D(:))*10)/10  ceil(max(D(:))*10)/10];
% % %   myfig('Communication Efficiency vs dist',ppostmp); B=bar(D',0.8); hold on; font(fntsz,usefont);
% % %   ylabel('Cohen''s d'); grid on; set(gca,'XGrid','off'); ylim(YLIM); 
% % %   B(1).FaceColor=[0 0 0]; B(1).FaceAlpha=1; B(2).FaceColor=[0 0 0]; B(2).FaceAlpha=.3;
% % %   for ii=1:Ncm-1; line([ii+.5 ii+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % %   set(gca,'XTickLabels',CMLBL_simple); 


% ==========  BOX + BEESWARM binned-scatter: CALIBER vs myelin: ED separated by 1st & 3rd quartiles 
% Switched order of caliber & myelin
  % ppostmp=getPlotPos([-2559 443 305 354],Ncm,0); ppostmp2=getPlotPos([-2559 9 305 354],Ncm,0); % big
  % ppostmp=getPlotPos([-2559 496 272 301],Ncm,0); ppostmp2=getPlotPos([-2559 115 272 301],Ncm,0); % small
  ppostmp=getPlotPos([-2559 548 230 249],Ncm,0); ppostmp2=getPlotPos([-2559 185 230 249],Ncm,0); % smaller
  scatopt={'filled','MarkerFaceAlpha',0.8,'MarkerEdgeColor','w','LineWidth',0.1};
  TLTOPT={'FontWeight','normal','FontSize',20};
  D3=EDnrm; Q1=prctile(D3(D3~=0),25); Q3=prctile(D3(D3~=0),75);
  mask1=D3>Q3; mask2=D3<Q1; mask1=mask1(~maskmd(:)); mask2=mask2(~maskmd(:));
  tstats=zeros(2,Ncm); cohens_d=zeros(2,Ncm);


  swarmopt.nbins = 70;                                               % more bins = finer density detail
  swarmopt.maxHalfWidth = 0.32;                                      % max horizontal half-width of swarm at peak density
  swarmopt.maxPts = 2000;                                            % points overlaid per group (set 0 to disable)
  swarmopt.tailPct  = 0.1;                                           % force include 0.5% lowest & highest (≈200 each tail for 40k)
  swarmopt.tailMax  = 80;                                          % cap each tail
  swarmopt.midBias  = 0;                                             % keep middle random (or set 1–2 to favor extremes more)
  swarmopt.ptSize = 4;                                               % marker size for subsampled points
  swarmopt.lineWidth = 2.0;                                          % thickness of density "bars"
  swarmopt.alphaPts = 0.5;                                           % marker alpha
  swarmopt.doColorPts = true;                                        % color subsampled points by y (turn off to save file size)
  swarmopt.seed = 1;                                                 % reproducible subsampling    


  for xx=1:Ncm 
      Dt1=ALLCM{2}{xx}; Dt2=ALLCM{1}{xx};                                   % zscored data
      X=Dt1(~maskmd(:)); Y=Dt2(~maskmd(:));                                 % All edges excluding main diagonal

    % ----- ED > Q3
    % Prep
      D1=X(mask1);   D2=Y(mask1);
      idx=unique([find(isnan(D1))   find(isnan(D2))]);
      D1(idx)=[];          D2(idx)=[];
      DPLT=[D1; D2];                                                        % Combine data into one vector
      tlbls=[]; for xxx=1:2; tlbls=[tlbls; xxx*ones(length(D1),1)]; end       % Create group labels
    % t-test
      [h,p,ci,stats]=ttest2(D1,D2);
    % Cohen's d
      t_m1=mean(D1);   t_m2=mean(D2);
      t_s1=std(D1);    t_s2=std(D2);
      t_n1=length(D1); t_n2=length(D2);
      t_pooled_sd=sqrt(((t_n1-1)*t_s1^2+(t_n2-1)*t_s2^2)/(t_n1+t_n2-2));
      t_cohens_d=(t_m1-t_m2)/t_pooled_sd;
    % Plot
      myfig(['ED > Q3: ' CMLBL_simple{xx}],ppostmp(xx,:)); hold on; colormap(cm2); clim([-3 3]);
      yline([0 0],'Color',[0.5 0.5 0.5],'LineWidth',1);
    % Density-swarm (binned) + optional point subsample
      plot_binswarm_density(gca, 1, D1, swarmopt);                       % group 1: Myelin
      plot_binswarm_density(gca, 2, D2, swarmopt);                       % group 2: Caliber
    % Add boxplots on top (no symbols)
      B=boxplot(DPLT,tlbls,'Colors','k','Widths',0.4,'Symbol','');
      set(findobj(B,'Type','Line'),'LineWidth', 2); font(fntsz,usefont);
      set(gca,'XTick',1:2,'XTickLabel',{'Myelin' 'Caliber'},'XTickLabelRotation',0); 
      ylabel(CMLBL_simple{xx}); box on;
      grid on; set(gca,'XGrid','off'); set(gca,'TickLength',[0 0]);
      % title([sprintf('t=%.1f, d=%.2f',stats.tstat,t_cohens_d)],TLTOPT{:});
      title([sprintf('d = %.2f',t_cohens_d)],TLTOPT{:});      
    % Save stats
      tstats(1,xx)=stats.tstat; cohens_d(1,xx)=t_cohens_d;
    % ----- ED < Q1
    % Prep
      D1=X(mask2);   D2=Y(mask2);
      idx=unique([find(isnan(D1))   find(isnan(D2))]);
      D1(idx)=[];          D2(idx)=[];
      DPLT=[D1; D2];                                                        % Combine data into one vector
      tlbls=[]; for xxx=1:2; tlbls=[tlbls; xxx*ones(length(D1),1)]; end       % Create group labels
    % t-test
      [h,p,ci,stats]=ttest2(D1,D2);
    % Cohen's d
      t_m1=mean(D1);   t_m2=mean(D2);
      t_s1=std(D1);    t_s2=std(D2);
      t_n1=length(D1); t_n2=length(D2);
      t_pooled_sd=sqrt(((t_n1-1)*t_s1^2+(t_n2-1)*t_s2^2)/(t_n1+t_n2-2));
      t_cohens_d=(t_m1-t_m2)/t_pooled_sd;
    % Plot
      myfig(['ED < Q1: ' CMLBL_simple{xx}],ppostmp2(xx,:)); hold on; colormap(cm2); clim([-3 3]);
      yline([0 0],'Color',[0.5 0.5 0.5],'LineWidth',1);
    % Density-swarm (binned) + optional point subsample
      plot_binswarm_density(gca, 1, D1, swarmopt);                       % group 1: Myelin
      plot_binswarm_density(gca, 2, D2, swarmopt);                       % group 2: Caliber
    % Add boxplots on top (no symbols)
      B=boxplot(DPLT,tlbls,'Colors','k','Widths',0.4,'Symbol','');
      set(findobj(B,'Type','Line'),'LineWidth', 2); font(fntsz,usefont);
      set(gca,'XTick',1:2,'XTickLabel',{'Myelin' 'Caliber'},'XTickLabelRotation',0); 
      ylabel(CMLBL_simple{xx}); box on;
      grid on; set(gca,'XGrid','off'); set(gca,'TickLength',[0 0]);
      % title([sprintf('t=%.1f, d=%.2f',stats.tstat,t_cohens_d)],TLTOPT{:});
      title([sprintf('d = %.2f',t_cohens_d)],TLTOPT{:});      
    % Save stats
      tstats(2,xx)=stats.tstat; cohens_d(2,xx)=t_cohens_d;
  end
% Summary: t-stats
  D=tstats; ppostmp=[-2559 600 569 197]; % ppostmp=[-2559 537 687 260];
  YLIM=[floor(min(D(:)))  ceil(max(D(:)))]; % YLIM=[floor(min(D(:))*10)/10  ceil(max(D(:))*10)/10];
  myfig('Communication Efficiency vs dist',ppostmp); B=bar(D',0.8); hold on; font(fntsz,usefont);
  ylabel('t-stat'); grid on; set(gca,'XGrid','off'); ylim(YLIM); 
  B(1).FaceColor=[0 0 0]; B(1).FaceAlpha=1; B(2).FaceColor=[0 0 0]; B(2).FaceAlpha=.3;
  for ii=1:Ncm-1; line([ii+.5 ii+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
  set(gca,'XTickLabels',CMLBL_simple); 
% Summary: Cohens d
  D=cohens_d; ppostmp=[-2559 600 569 197]; %ppostmp=[-2559 537 687 260];
  YLIM=[-1 1]; %YLIM=[floor(min(D(:))*10)/10  ceil(max(D(:))*10)/10];
  myfig('Communication Efficiency vs dist',ppostmp); B=bar(D',0.8); hold on; font(fntsz,usefont);
  ylabel('Cohen''s d'); grid on; set(gca,'XGrid','off'); ylim(YLIM); 
  B(1).FaceColor=[0 0 0]; B(1).FaceAlpha=1; B(2).FaceColor=[0 0 0]; B(2).FaceAlpha=.3;
  for ii=1:Ncm-1; line([ii+.5 ii+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
  set(gca,'XTickLabels',CMLBL_simple); set(gca,'TickLength',[0 0]);

% TMP SAVE
  pp=1; set(gcf,'Renderer','painters'); svstr=[num2str(pp) '_gt-Q3_' CMLBL_simple{pp}];
  exportgraphics(gcf,['~/Downloads/3-' svstr  '_nolbl.pdf'],'ContentType','vector');




% % % % ========== BARPLOT: SA-axis vs Myelin_CMi - Caliber_CMi
% % % % --- Compute correlation of ΔCommunication with SAaxis
% % %   svcorr=zeros(1,Ncm); pvals=zeros(1,Ncm);
% % %   for xx=1:Ncm 
% % %       Dt1=ALLCM{2}{xx}; Dt2=ALLCM{1}{xx};                                % zscored data
% % %       MV=abs(min(Dt1(:))); MV=MV+(.001*MV); Dt1=Dt1+MV;                  % Shift all myelin values to > 0
% % %       MV=abs(min(Dt2(:))); MV=MV+(.001*MV); Dt2=Dt2+MV;                  % Shift all caliber values to > 0
% % %       Dt1=nanmean(Dt1); Dt2=nanmean(Dt2);                                % Compute mean communication for each node
% % %       Dx=Dt1-Dt2;                                                        % Difference of myelin & caliber
% % %       [svcorr(xx),pvals(xx)]=corr(Dx',SA_map);                           % Correlate with SAaxis
% % %   end
% % % % --- Plot
% % %   D=svcorr; ppostmp=[-2559 564 357 233]; YLIM=[floor(min(D(:))*10)/10  ceil(max(D(:))*10)/10];
% % %   myfig('ΔCommunication vs SA-axis',ppostmp); B=bar(D,0.8); hold on; font(fntsz,usefont);
% % %   grid on; set(gca,'XGrid','off'); ylim(YLIM); % ylabel('r(ΔCommunication, SA-axis)');
% % %   B.CData=cmap_6cm; B.FaceColor='flat';
% % %   set(gca,'XTickLabels',CMLBL_simple);
% % % % overlay significance markers
% % %   for xx=1:Ncm
% % %         if pvals(xx)<0.05
% % %             text(xx,svcorr(xx)+0.03,'*','HorizontalAlignment','center','FontSize',30);
% % %         end
% % %   end
% % % 
% % % 
% % % %% ----- COMPARE TO FC ACROSS DISTANCE BINS
% % % 
% % % % ========== Matrices: corr(myelin_CMi, FCx) across distance bins
% % % % Create ED mask
% % %   D=EDnrm; Q1=prctile(D(D~=0),25); Q2=prctile(D(D~=0),50); Q3=prctile(D(D~=0),75);
% % %   mask=double(D<Q1); mask(D>Q1 & D<Q2)=2; mask(D>Q2 & D<Q3)=3; mask(D>Q3)=4; mask(logical(eye(Nnode)))=0;
% % % % Correlated across distance bins
% % %   svcorr=zeros(Ncm,4,Nfc); pvals=zeros(Ncm,4,Nfc);
% % %   for xx=1:Ncm
% % %       Dt1=ALLCM{2}{xx};                                                    % z-scored myelin communication data
% % %       for yy=1:Nfc
% % %             Dt2=Yz_m{yy};
% % %           % Across distance bins
% % %             for zz=1:4
% % %                 Dx=Dt1(mask==zz); Dy=Dt2(mask==zz);
% % %                 Dy(isnan(Dx))=[]; Dx(isnan(Dx))=[];                      % remove null values
% % %                 [svcorr(xx,zz,yy),pvals(xx,zz,yy)]=corr(Dx,Dy,'type','Spearman');
% % %             end
% % %       end
% % %   end
% % % % Plot: matrices
% % %   ppostmp=getPlotPos([-2558 495 252 338],Nfc,0);
% % %   CMAX=ceil(max(abs(svcorr(:)))*10)/10; CLIM=[-1*CMAX  CMAX];
% % %   for yy=1:Nfc
% % %        myfig(['Rho(Myelin-CM, ' ylbl{yy} ')'],ppostmp(yy,:)); 
% % %        imagesc(svcorr(:,:,yy)); colormap(cm2); clim(CLIM); % title(ylbl{yy});
% % %        set(gca,'YTickLabels',CMLBL_simple,'XTickLabels',''); xlabel('Distance bins'); font(fntsz,usefont);
% % %   end
% % % % Plot: lines
% % %   ppostmp=getPlotPos([-2558 41 252 338],Nfc,0);
% % %   YLIM=[floor(min(svcorr(:))*10)/10  ceil(max(svcorr(:))*10)/10];
% % %   for yy=1:Nfc
% % %       myfig(['Rho(Myelin-CM, ' ylbl{yy} ')'],ppostmp(yy,:)); hold on;
% % %       for xx=1:Ncm
% % %             plot(1:4,squeeze(svcorr(xx,:,yy)),'--','LineWidth',3,'Color',cmap_6cm(xx,:));
% % %       end
% % %       xlabel('Distance bins'); font(fntsz,usefont); % ylabel('Spearman''s \rho');
% % %       set(gca,'XTickLabels',''); grid on; box on; ylim(YLIM); % title(ylbl{yy});
% % %       line(xlim,[0 0],'Color','k','LineWidth',1);
% % %   end
% % % 
% % % 
% % % % ========== Regression modeling of FC with myelin_CM + interaction with ED
% % %   Xd=EDnrm; Xdz=nzzscore(Xd,0); 
% % %   ylbl_simple=strrep(ylbl,'_',''); mdlcoeffs=zeros(Ncm,Nfc,5); pvals=zeros(Ncm,Nfc,5);
% % %   for xx=1:Ncm
% % %       Dt2=ALLCM{2}{xx}; Dt1=ALLCM{1}{xx};                                 % z-scored myelin & caliber communication data
% % %       for yy=1:Nfc
% % %             Dt3=Yz_m{yy};
% % %           % Prep data
% % %             tbl = table(Dt3(:), Dt2(:), Dt1(:), Xdz(:),'VariableNames',{ylbl_simple{yy},'myelin','caliber','ED'});
% % %             tbl = tbl(~any(ismissing(tbl) | tbl{:,:} == 0, 2), :);        % Remove rows containing NaNs or 0s in any variable
% % %             tbl.Int = tbl.myelin .* tbl.ED;
% % %           % Model
% % %             mdl = fitlm(tbl, [ylbl_simple{yy} ' ~ myelin + caliber + ED + Int']);
% % %           % Extract coefficients
% % %             mdlcoeffs(xx,yy,:) = mdl.Coefficients.Estimate;
% % %             pvals(xx,yy,:)     = mdl.Coefficients.pValue;
% % %       end
% % %   end
% % % % --- Plot: Matrix
% % %   D=mdlcoeffs(:,:,5); D2=pvals(:,:,5); D(D2>0.05)=0;
% % %   CMAX=ceil(max(abs(D(:)))*10)/10; CLIM=[-1*CMAX  CMAX];
% % %   myfig('MyelinCMi*ED Betas',[-2558 410 477 387]); imagesc(D); clim(CLIM); colormap(cm2);
% % %   set(gca,'YTickLabels',CMLBL_simple,'XTick',[1:Nfc],'XTickLabels',ylbl_simple); font(fntsz,usefont);
% % % % --- Plot: bars
% % %   D=mdlcoeffs(:,:,5); D2=pvals(:,:,5); D(D2>0.05)=0;
% % %   YLIM=[floor(min(D(:))*10)/10  ceil(max(D(:))*10)/10];
% % %   myfig('MyelinCMi*ED Betas',[-2046 410 922 387]); B=bar(D',0.8);
% % %   ylabel('\beta myelinCM*ED','FontWeight','bold');
% % %   for xx=1:Ncm; B(xx).FaceColor=cmap_6cm(xx,:); end
% % %   grid on; set(gca,'XGrid','off'); ylim(YLIM);
% % %   for ii=1:Nfc-1; line([ii+.5 ii+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % %   set(gca,'XTickLabels',ylbl_simple); font(fntsz,usefont);








% -------------------------------------------------------------------------
%% ----- COMPARE SPECTRAL PROPERTIES: vs Caliber
% Extend to include g-ratio
  kmax=50; doSignAlign=true; 
  D1=SC; [Ahat_D1,Lsym_D1] = normalize_operators(D1);

% MTsat VS Caliber
  D2=MC; tmpstr='caliber vs MTsat';
  [Ahat_D2,Lsym_D2] = normalize_operators(D2);
  spec_out_myelin=spectral_alignment_compare_adj_lap(Ahat_D1,Ahat_D2,kmax,doSignAlign,tmpstr);

% g-ratio VS Caliber
  D2=MC2; tmpstr='caliber vs g-ratio';
  [Ahat_D2,Lsym_D2] = normalize_operators(D2);
  spec_out_gratio=spectral_alignment_compare_adj_lap(Ahat_D1,Ahat_D2,kmax,doSignAlign,tmpstr);

% Delay VS Caliber
  D2=r_delay; tmpstr='caliber vs Delay (rate)';
  [Ahat_D2,Lsym_D2] = normalize_operators(D2);
  spec_out_delay=spectral_alignment_compare_adj_lap(Ahat_D1,Ahat_D2,kmax,doSignAlign,tmpstr);

% --- Plot mean cosine similarity for myelin & delay
  myfig(['Cumulative spectral alignment with caliber: ' str],[-2559 439 915 358]);
% Adjacency
  subplot(1,2,1); 
  plot(1:kmax,spec_out_myelin.cos_sim_adj,'-','LineWidth',3,'Color',cmap_pred(2,:)); hold on;
  plot(1:kmax,spec_out_gratio.cos_sim_adj,'-','LineWidth',3,'Color',cmap_pred(3,:));
  plot(1:kmax,spec_out_delay.cos_sim_adj,'-','LineWidth',3,'Color',cmap_pred(4,:));
  xlabel('k (leading eigenvectors)'); ylabel('Mean cosine similarity');
  title('Adjacency'); ylim([0 1]); xlim([0 kmax]); grid on; font(fntsz,usefont); 
% Laplacian
  subplot(1,2,2); 
  plot(1:kmax,spec_out_myelin.cos_sim_lap,'-','LineWidth',3,'Color',cmap_pred(2,:)); hold on;
  plot(1:kmax,spec_out_gratio.cos_sim_lap,'-','LineWidth',3,'Color',cmap_pred(3,:));
  plot(1:kmax,spec_out_delay.cos_sim_lap,'-','LineWidth',3,'Color',cmap_pred(4,:));
  xlabel('k (leading eigenvectors)'); ylabel('Mean cosine similarity');
  title('Laplacian'); ylim([0 1]); xlim([0 kmax]); grid on; font(fntsz,usefont);

% Joint plot 1: 
% MTsat vs delay
  plot_spectral_summary_2(spec_out_myelin,spec_out_delay,cmap_pred([2 4],:),{'MTsat' 'delay'},'vs Caliber')
% MTsat vs g-ratio
  plot_spectral_summary_2(spec_out_myelin,spec_out_gratio,cmap_pred([2 3],:),{'MTsat' 'g-ratio'},'vs Caliber')
% All 3: MTsat, g-ratio, delay
  plot_spectral_summary_3(spec_out_myelin,spec_out_gratio,spec_out_delay,cmap_pred(2:4,:),{'MTsat' 'g-ratio' 'delay'},'vs caliber',[1 0 1])
  title ''; legend off; set(gca,'TickLength',[0 0]); subplot(1,2,1); legend off; title ''; set(gcf,'Position',[-2559 482 973 315]);
% Save
  svstr='1_alignment-distance';
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '_k50.pdf'],'ContentType','vector');


% Joint plot 2: Radar-style
% MTsat & delay
  plot_spectral_summary_radar(spec_out_myelin,spec_out_delay,cmap_pred([2 4],:),{'MTsat' 'Delay'},'vs Caliber')
% MTsat & g-ratio
  plot_spectral_summary_radar(spec_out_myelin,spec_out_gratio,cmap_pred([2 3],:),{'MTsat' 'g-ratio'},'vs Caliber')



%% ======== SPECTRAL ANALYSIS (WINDOWS)


% ======= WINDOWED SUBSPACES PART I
% Manual windows to get started

% Example windows from global→local scales
  kmax=399;
  wins = {2:6, 5:15, 14:40, 40:kmax};

% Row-normalized weights matrix for Morans I
  W_spatial = rowStandardizeWeights(BIN);

% Options
opts = struct(                  ...
  'kmax', kmax,                 ...
  'normalize', true,            ...
  'wins', {wins},               ...
  'doSignAlign', true,          ...
  'coords', pinfo.coor,         ...
  'SA_axis', SA_map,            ...
  'modules', YeoRSNs,           ...
  'Wspatial', W_spatial,        ...
  'plotMaps', false,            ...
  'cmap', cmap_pred(2:4,:),     ... 
  'matchMode', 'global-banded', ...
  'band', 3 );

% === Compare all myelin metrics to caliber
  tmpstr='caliber vs myelin metrics';
  D=zeros(Nnode,Nnode,3); D(:,:,1)=MC; D(:,:,2)=MC2; D(:,:,3)=r_delay;
  tmplbls={'MTsat' 'g-ratio' 'delay'};
  out = spectral_alignment_windows_j(SC,D,tmplbls,tmpstr,opts);             % Original
  set(out.figs.subspace,'Position',[-2559 476 664 321]); 
  % set(out.figs.subspace,'Position',[-2559 518 605 279]); % smaller
  set(out.figs.per_k,'Position',[-2559 -539 2560 520])

% === Compare topology of worst aligned eigs
  % idx_eigs_adj=[4 5 7 8 9 12 14 17 18 19 22 23 26 30 40 45 50];
  idx_eigs_adj=[4 5 7 8 9 12 14 17 18 19 22 23 26 30 40 45 50 56 61 85 86 92 98 109 134 146 161 166 176 184 192];
  plot_matched_eig_metrics(out,'adj',idx_eigs_adj,false,2,'both','norm','off'); % set(gcf,'Position',[-2559 497 1077 300]); % for bar plot
  plot_matched_eig_metrics(out,'lap',idx_eigs_adj,false,2,'both','norm','off'); % set(gcf,'Position',[-2559 497 1077 300]); % for bar plot



% ======= WINDOWED SUBSPACES PART II
% Auto-define windows


% ------- SPECTRAL FINGERPRINTS
  kmax=399;
  wins='auto';  
  % wins = {2:6, 7:30, 31:80, 81:200, 201:kmax};
% --- Build Lap basis separate to ensure done correctly
  W = double(SCnrm); W = (W+W.')/2;                                         % the reference structural matrix used to compute DE (e.g., caliber)
  s = sum(W,2); s = max(s, eps);                      % avoid zeros
  Dinvsqrt = spdiags(s.^(-0.5), 0, Nnode, Nnode);
  Dhalf    = spdiags(s.^( 0.5), 0, Nnode, Nnode);
  Lsym = speye(Nnode) - Dinvsqrt * W * Dinvsqrt;
  % % % D   = spdiags(s,0,Nnode,Nnode);
  % % % Lsym = speye(Nnode) - (D^(-1/2)) * W * (D^(-1/2));                  % Normalized Laplacian (symmetric)
  [Vlap, ~] = eigs(Lsym, kmax, 'smallestabs');                              % Leading K Laplacian eigenvectors (smallest eigenvalues)
  Vbasis.lap.V = Vlap;
  Vbasis.adj.V = out.adj.Vref(:,1:kmax);                                    % using already-computed ADJ basis
  % Vbasis.lap.V = out.lap.Vref(:,1:kmax);
  Vbasis.lap.V = Vlap;
  OpForModel = {'adj','adj','adj','adj','adj','lap'};                       % Operator mapping per model:
  dsLabels = {'caliber','MTsat','g-ratio','delay'};                         % Labels:
  useCM=[ALLCMnx {ALLCMnx_delay}];  useCMLBL=CMLBL_simple;                  % Combine ALLCM & LBLS
% Run:
  res=spectral_fingerprint_windows(useCM,useCMLBL,Vbasis,wins, ...
      'Lsym', Lsym,                       ...
      'QuantileEdges', [0 .05 .25 .55 1], ...
      'StartAt', 2,                       ...
      'RefIdx',1,                         ... 
      'OpForModel',OpForModel,            ...
      'DatasetLabels',dsLabels,           ... 
      'Cmap', cmap_pred(1:4,:),           ...
      'LabelMode','off',                   ...
      'DE_UseDhalf',true,                 ...
      'Dhalf',Dhalf);

% set(gcf,'Position',[-2559 460 1508 337]) % main plot
% set(gcf,'Position',[-1050 625 781  172]) % k=1 plot

% TMP SAVE
  set(gcf,'Renderer','painters'); svstr='1_fingerprints-k1';
  exportgraphics(gcf,['~/Downloads/' svstr  '.pdf'],'ContentType','vector');

  set(gcf,'Renderer','painters'); svstr='2_fingerprints-k2-399';
  exportgraphics(gcf,['~/Downloads/' svstr  '.pdf'],'ContentType','vector');



% ------- WINDOW SUMMARIES & MATCHED EIGS
% Using windows: rerun summary
  wins=res.wins;              % use windows from fingerprints function

% Options
opts = struct(                  ...
  'kmax', kmax,                 ...
  'normalize', true,            ...
  'wins', {wins},               ...
  'doSignAlign', true,          ...
  'coords', pinfo.coor,         ...
  'SA_axis', SA_map,            ...
  'modules', YeoRSNs,           ...
  'Wspatial', W_spatial,        ...
  'plotMaps', false,            ...
  'cmap', cmap_pred(2:4,:),     ... 
  'matchMode', 'global-banded', ...
  'LabelMode', 'off',            ...
  'KPerPage', 100,              ...
  'band', 3 );

% === Compare all myelin metrics to caliber: ROBUST version
  tmpstr='caliber vs myelin metrics';
  D=zeros(Nnode,Nnode,3); D(:,:,1)=MC; D(:,:,2)=MC2; D(:,:,3)=r_delay;
  tmplbls={'MTsat' 'g-ratio' 'delay'};
  % out = spectral_alignment_windows_j(SC,D,tmplbls,tmpstr,opts);
  out = spectral_alignment_windows_robust(SC,D,tmplbls,tmpstr,opts);        % ROBUST 
  set(out.figs.subspace,'Position',[-2559 550 622 247]);    % set(out.figs.subspace,'Position',[-2559 476 664 321]);  
  set(out.figs.per_k,'Position',[-2559 -539 2560 520])
  % AX=findall(gcf,'type','axes'); set(AX(1),'TickLength',[0 0]);
% Save
  svstr='2_window-alignment';
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '.pdf'],'ContentType','vector');


% === Compare topology of worst aligned eigs
  % idx_eigs_adj=[7 8 9 12 14 17 18 19 22 23 26 30 40 45 50];
  idx_eigs_adj=[7 8 9 12 14 17 18 19 22 23 26 30 40 45 50 56 61 85 86 92 98 ... 
                109 134 146 161 166 176 184 192 ... 
                201 213 228 248 249 259 261 272 275 280 281 ...
                313 314 316 317 319 337 341 349 350 358 368 378 379 383 385 395 396 397 399];
  plot_matched_eig_metrics(out,'adj',idx_eigs_adj,false,2,'bar','norm','off','metric'); % set(gcf,'Position',[-2559 497 1157 300]); % for bar plot
  plot_matched_eig_metrics(out,'lap',idx_eigs_adj,false,2,'bar','norm','off','metric'); % set(gcf,'Position',[-2559 497 1157 300]); % for bar plot

% Save
  % svstr='3_matched-eigs-bymetric_1_adj';
  svstr='3_matched-eigs-bymetric_2_lap';
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '.pdf'],'ContentType','vector');



% === Summarize topology across all k
% --- Boxplot
% % % metrics = {'delta_SAr','delta_smooth','delta_PR','delta_modR2','delta_MoranI'};
% % % figs_adj = delta_topology_boxplots(out,'adj',metrics,tmplbls,out.cmap);
% % % figs_lap = delta_topology_boxplots(out,'lap',metrics,tmplbls,out.cmap);
% --- Boxchart
  metrics = {'delta_SAr','delta_smooth','delta_PR','delta_modR2','delta_MoranI'};
  ppostmp=getPlotPos([-2559 546 648 251],numel(metrics),0);
% Adj
  figs_adj = delta_topology_boxcharts(out,'adj',metrics,tmplbls,out.cmap,'off'); 
  for ff=1:numel(metrics); set(figs_adj(ff),'Position',ppostmp(1,:)); end   % overlapping
  % for ff=1:numel(metrics); set(figs_adj(ff),'Position',ppostmp(ff,:)); end  % non-overlapping
% Lap
  figs_lap = delta_topology_boxcharts(out,'lap',metrics,tmplbls,out.cmap,'off');
  for ff=1:numel(metrics); set(figs_lap(ff),'Position',ppostmp(1,:)); end   % overlapping
  % for ff=1:numel(metrics); set(figs_lap(ff),'Position',ppostmp(ff,:)); end  % non-overlapping






% --- Check range of kmax
% % %   Ks = [120 200 300 399];
% % %   all_res = cell(size(Ks));
% % % % - Build Lap basis separate to ensure done correctly
% % %   W = double(SCnrm); W = (W+W.')/2;                                           % the reference structural matrix used to compute DE (e.g., caliber)
% % %   s = sum(W,2); s = max(s, eps);                                           % avoid zeros
% % %   Dinvsqrt = spdiags(s.^(-0.5), 0, Nnode, Nnode);
% % %   Dhalf    = spdiags(s.^( 0.5), 0, Nnode, Nnode);
% % %   Lsym = speye(Nnode) - Dinvsqrt * W * Dinvsqrt;
% % %   [Vlap, ~] = eigs(Lsym, kmax, 'smallestabs');                                % Leading K Laplacian eigenvectors (smallest eigenvalues)
% % % % Operator mapping per model:
% % %   OpForModel = {'adj','adj','adj','adj','adj','lap'};                       % last is DE
% % % % Combine ALLCM & LBLS
% % %   useCM=[ALLCMnx {ALLCMnx_delay}];  useCMLBL=CMLBL_simple;  %useCMLBL=[CMLBLnx {CMLBLnx_delay}] ;
% % % % loop
% % %   for t=1:numel(Ks)
% % %       K = Ks(t);
% % %       Vbasis.adj.V = out.adj.Vref(:,1:K);
% % %       Vbasis.lap.V = Vlap(:,1:K);                                           % your Lsym basis
% % %       resK = spectral_fingerprint_windows(useCM, useCMLBL, Vbasis, wins_for(K), ...
% % %             'RefIdx',1,'OpForModel',OpForModel,'Dhalf',Dhalf,'DE_UseDhalf',true,'Plot','none');
% % %       all_res{t} = resK;
% % %   end
% % % % Compare to the largest-K as reference:
% % %   ref = all_res{end}.phi;
% % %   for t=1:numel(Ks)-1
% % %       diffmax(t) = max(abs(all_res{t}.phi(:) - ref(:)), [], 'omitnan');
% % %   end
% % %   table(Ks(1:end-1).', diffmax.', 'VariableNames',{'K','maxAbsDiff'})



% -------------------------------------------------------------------------






%% ----- COMPARE SPECTRAL PROPERTIES: vs caliber
  kmax=50; doSignAlign=true; 
  D1=SC; [Ahat_D1,Lsym_D1] = normalize_operators(D1);

% MTsat VS Caliber
  D2=MC; tmpstr='caliber vs MTsat';
  [Ahat_D2,Lsym_D2] = normalize_operators(D2);
  spec_out_myelin=spectral_alignment_compare_adj_lap(Ahat_D1,Ahat_D2,kmax,doSignAlign,tmpstr);

% Delay VS Caliber
  D2=r_delay; tmpstr='caliber vs Delay (rate)';
  [Ahat_D2,Lsym_D2] = normalize_operators(D2);
  spec_out_delay=spectral_alignment_compare_adj_lap(Ahat_D1,Ahat_D2,kmax,doSignAlign,tmpstr);

% --- Plot mean cosine similarity for myelin & delay
  myfig(['Cumulative spectral alignment with caliber: ' str],[-2559 439 915 358]);
% MTsat
  subplot(1,2,1); 
  plot(1:kmax,spec_out_myelin.cos_sim_adj,'-','LineWidth',3,'Color',cmap_pred(2,:)); hold on;
  plot(1:kmax,spec_out_delay.cos_sim_adj,'-','LineWidth',3,'Color',cmap_pred(4,:));
  xlabel('k (leading eigenvectors)'); ylabel('Mean cosine similarity');
  title('Adjacency'); ylim([0 1]); xlim([0 kmax]); grid on; font(20,'Cambria'); 
% delay
  subplot(1,2,2); 
  plot(1:kmax,spec_out_myelin.cos_sim_lap,'-','LineWidth',3,'Color',cmap_pred(2,:)); hold on;
  plot(1:kmax,spec_out_delay.cos_sim_lap,'-','LineWidth',3,'Color',cmap_pred(4,:));
  xlabel('k (leading eigenvectors)'); ylabel('Mean cosine similarity');
  title('Laplacian'); ylim([0 1]); xlim([0 kmax]); grid on; font(20,'Cambria');

% Joint plot 1
  plot_spectral_summary(spec_out_myelin,spec_out_delay,cmap_pred([2 4],:),{'MTsat' 'Delay'},'vs Caliber')
% Joint plot 2: Radar-style
  plot_spectral_summary_radar(spec_out_myelin,spec_out_delay,cmap_pred([2 4],:),{'MTsat' 'Delay'},'vs Caliber')


% ======= ADD-ONs: 1 (individual targets)
% windowed-subspace analysis of Eigs between Adj & Lap
% matching of Eigs within windows & topology comparison
% kernel-weighted contributions of windows

% Example windows from global→local scales
  wins = {2:3, 4:8, 9:12, 13:16, 17:35, 36:50};

% Row-normalized weights matrix for Morans I
  W_spatial = rowStandardizeWeights(BIN);

% Options
opts = struct( ...
  'kmax', 50, ...
  'normalize', true, ...
  'wins', {wins}, ...
  'doSignAlign', true, ...
  'coords', pinfo.coor, ...
  'SA_axis', SA_map, ...
  'modules', YeoRSNs, ...
  'Wspatial', W_spatial, ...
  'plotMaps', false, ...
  'kernel', struct('doAdj', true, 'doLap', true, 'beta', 1, 'tau', 1) );

% Compare MTsat vs CALIBER
  tmpstr='caliber vs MTsat';
  out_mtsat = analyze_spectral_alignment_windows(SC,MC,tmpstr,opts);

% Compare DELAY vs CALIBER
  tmpstr='caliber vs DELAY (rate)';
  out_delay = analyze_spectral_alignment_windows(SC,r_delay,tmpstr,opts);





% -------------------------------------------------------------------------



% % % % ==========   Scatter plot caliber vs myelin: data colored by ED
% % %   xx=6; 
% % %   ppostmp=[-2558 437 415 360];
% % % % ---- zscored data
% % %   D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; % caliber communication model
% % %   D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx}; % myelin communication model
% % %   D3=EDnrm;
% % % % -- All edges colored by connectivity
% % %   TMPLBL2='caliber vs myelin (data colored by Euclid dist)';
% % %   X=D1(:); Y=D2(:); Z=D3(:);                                                % flatten
% % %   X(maskmd(:))=[]; Y(maskmd(:))=[]; Z(maskmd(:))=[];                        % Remove main diagonal
% % % % Plot
% % %   myfig(TMPLBL2,ppostmp); hold on;
% % %   scatter(X,Y,5,Z,'filled','MarkerFaceAlpha', 0.7);
% % %   colormap(cm3); CB=colorbar; %ylabel(CB,'ED','FontSize',fntsz);
% % %   xlabel(extractBefore(Dstr1,'-')); ylabel(extractBefore(Dstr2,'-')); title(CMLBL_simple{xx}); font(fntsz,usefont);
% % %   box on;
% % %   % ylim([-3 5]); xlim([-2.1 6]); % SPE
% % %   % xlim([-3.9 6]); ylim([-5 5]); % NE
% % %   % ylim([-1.6 3.6]); xlim([-2.2 6]); % SIE
% % %   % xlim([-2.5 4.2]); ylim([-2.5 4.5]); % PT
% % %   % xlim([-3.5 5]); ylim([-3 4]); % DE
% % %   line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','--','LineWidth',3);
% % %   line([0 0],ylim,'Color',[.5 .5 .5],'LineStyle','--','LineWidth',3);
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % ==========   Visualize CMs general trends
% % % % Caliber vs myelin for each CM
% % % % vs BOLD-FC
% % % % vs ED
% % % 
% % %   xx=1; 
% % %   opt={'XTickLabels','','YTickLabels',''}; opt4={'XTickLabel',{'all' 'dis' 'con'},'XTickLabelRotation',45};
% % %   opt2={'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk}; opt3={'type','Spearman'}; 
% % %   cmaptmp=cmap_6cm(xx,:); 
% % %   ppostmp=getPlotPos([-2559 300 564 497],4);   ppostmp2=getPlotPos([-2559 -13 363 309],6); 
% % %   ppostmp3=getPlotPos([-2559 -275 333 261],7); ppostmp4=getPlotPos([-2559 -539 333 261],7);
% % % 
% % % % ---- zscored data
% % %   D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; CMAX=prctile(abs(D1(D1~=0)),99); CLIM1=[CMAX*-1 CMAX]; %CMAX=abs(max(D1(:)));
% % %   D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx}; CMAX=prctile(abs(D2(D2~=0)),99); CLIM2=[CMAX*-1 CMAX]; %CMAX=abs(max(D2(:)));
% % % % -- Correlations with ED
% % %   D3=ED; YSTR='ED'; svcorr=zeros(1,6); useppos=ppostmp4; TMPLBL2=[CMLBL_simple{xx} ' vs ' YSTR];
% % %   X=D1;          Y=D3;          ii=1; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (ALL)'],[YSTR ' (ALL)'],str,cm,useppos(ii,:),20); % D1 all
% % %   X=D2;          Y=D3;          ii=2; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (ALL)'],[YSTR ' (ALL)'],str,cm,useppos(ii,:),20); % D2 all
% % %   X=D1(~masknz); Y=D3(~masknz); ii=3; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (DIS)'],[YSTR ' (DIS)'],str,cm,useppos(ii,:),20); % D1 DIS
% % %   X=D2(~masknz); Y=D3(~masknz); ii=4; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (DIS)'],[YSTR ' (DIS)'],str,cm,useppos(ii,:),20); % D2 DIS
% % %   X=D1(masknz);  Y=D3(masknz);  ii=5; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (CON)'],[YSTR ' (CON)'],str,cm,useppos(ii,:),20); % D1 CONN
% % %   X=D2(masknz);  Y=D3(masknz);  ii=6; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (CON)'],[YSTR ' (CON)'],str,cm,useppos(ii,:),20); % D2 CONN
% % % % Summary bar plot
% % %   myfig([YSTR str],useppos(7,:)); B=bar(svcorr,0.6);  hold on; title(TMPLBL2);
% % %   ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %   for pp=1:2:6; B.CData(pp,:)=cmap_pred(1,:); end 
% % %   for pp=2:2:6; B.CData(pp,:)=cmap_pred(2,:); end
% % %   B.FaceColor='flat'; set(gca,'XTick',1.5:2:5.5,opt4{:}); 
% % %   for pp=2:2:4; line([pp+.5 pp+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % % % -- Correlations with BOLD-FC
% % %   D3=FC; YSTR='BOLD-FC'; svcorr=zeros(1,6); useppos=ppostmp3; TMPLBL2=[CMLBL_simple{xx} ' vs ' YSTR];
% % %   X=D1;          Y=D3;          ii=1; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (ALL)'],[YSTR ' (ALL)'],str,cm,useppos(ii,:),20); % D1 all
% % %   X=D2;          Y=D3;          ii=2; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (ALL)'],[YSTR ' (ALL)'],str,cm,useppos(ii,:),20); % D2 all
% % %   X=D1(~masknz); Y=D3(~masknz); ii=3; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (DIS)'],[YSTR ' (DIS)'],str,cm,useppos(ii,:),20); % D1 DIS
% % %   X=D2(~masknz); Y=D3(~masknz); ii=4; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (DIS)'],[YSTR ' (DIS)'],str,cm,useppos(ii,:),20); % D2 DIS
% % %   X=D1(masknz);  Y=D3(masknz);  ii=5; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (CON)'],[YSTR ' (CON)'],str,cm,useppos(ii,:),20); % D1 CONN
% % %   X=D2(masknz);  Y=D3(masknz);  ii=6; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (CON)'],[YSTR ' (CON)'],str,cm,useppos(ii,:),20); % D2 CONN
% % % % Summary bar plot
% % %   myfig([YSTR str],useppos(7,:)); B=bar(svcorr,0.6);  hold on; title(TMPLBL2);
% % %   ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %   for pp=1:2:6; B.CData(pp,:)=cmap_pred(1,:); end 
% % %   for pp=2:2:6; B.CData(pp,:)=cmap_pred(2,:); end
% % %   B.FaceColor='flat'; set(gca,'XTick',1.5:2:5.5,opt4{:}); 
% % %   for pp=2:2:4; line([pp+.5 pp+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % % % Histograms
% % %   myfig(str,ppostmp2(1,:)); histogram(D1,bins); title(Dstr1); font(fntsz,usefont);
% % %   myfig(str,ppostmp2(2,:)); histogram(D2,bins); title(Dstr2); font(fntsz,usefont);
% % % % -- Correlations with each other
% % %   svcorr=zeros(1,3); useppos=ppostmp2; TMPLBL2=['caliber vs myelin (' CMLBL_simple{xx} ')'];
% % %   X=D1;          Y=D2;          [~,svcorr(1)]=ddisp_corr(X,Y,[Dstr1 ' (ALL)'],[Dstr2 ' (ALL)'],str,cm,useppos(3,:),20);
% % %   X=D1(~masknz); Y=D2(~masknz); [~,svcorr(2)]=ddisp_corr(X,Y,[Dstr1 ' (DIS)'],[Dstr2 ' (DIS)'],str,cm,useppos(4,:),20);
% % %   X=D1(masknz);  Y=D2(masknz);  [~,svcorr(3)]=ddisp_corr(X,Y,[Dstr1 ' (CON)'],[Dstr2 ' (CON)'],str,cm,useppos(5,:),20);
% % % % Summary bar plot
% % %   myfig([TMPLBL2 str],useppos(6,:)); B=bar(svcorr,0.6);  hold on; title(TMPLBL2);
% % %   ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %   B.CData=cmaptmp; B.FaceColor='flat'; set(gca,opt4{:}); 
% % % % Connectomes: sequential node ordering
% % %   connplot(D1,[TMPLBL2 str],[],ppostmp(1,:)); title(Dstr1); clim(CLIM1); colormap(cm2); set(gca,opt{:})
% % %   connplot(D2,[TMPLBL2 str],[],ppostmp(2,:)); title(Dstr2); clim(CLIM2); colormap(cm2); set(gca,opt{:})
% % % % Connectomes: Yeo RSN ordering
% % %   connplot(D1,[TMPLBL2 str],pinfo,ppostmp(3,:)); title(Dstr1); clim(CLIM1); colormap(cm2); set(gca,opt2{:})
% % %   connplot(D2,[TMPLBL2 str],pinfo,ppostmp(4,:)); title(Dstr2); clim(CLIM2); colormap(cm2); set(gca,opt2{:})
% % % 
% % % 
% % % % ---- NO xfm data
% % %   D1=ALLCMnx{1}{xx}; Dstr1=CMLBLnx{1}{xx}; CLIM1=[prctile(D1(D1~=0),1)  prctile(D1(D1~=0),99)];
% % %   D2=ALLCMnx{2}{xx}; Dstr2=CMLBLnx{2}{xx}; CLIM2=[prctile(D2(D2~=0),1)  prctile(D2(D2~=0),99)];
% % % % -- Correlations with ED
% % %   D3=ED; YSTR='ED'; svcorr=zeros(1,6); useppos=ppostmp4; TMPLBL2=[CMLBL_simple{xx} ' vs ' YSTR];
% % %   X=D1;          Y=D3;          ii=1; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (ALL)'],[YSTR ' (ALL)'],str,cm,useppos(ii,:),20); % D1 all
% % %   X=D2;          Y=D3;          ii=2; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (ALL)'],[YSTR ' (ALL)'],str,cm,useppos(ii,:),20); % D2 all
% % %   X=D1(~masknz); Y=D3(~masknz); ii=3; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (DIS)'],[YSTR ' (DIS)'],str,cm,useppos(ii,:),20); % D1 DIS
% % %   X=D2(~masknz); Y=D3(~masknz); ii=4; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (DIS)'],[YSTR ' (DIS)'],str,cm,useppos(ii,:),20); % D2 DIS
% % %   X=D1(masknz);  Y=D3(masknz);  ii=5; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (CON)'],[YSTR ' (CON)'],str,cm,useppos(ii,:),20); % D1 CONN
% % %   X=D2(masknz);  Y=D3(masknz);  ii=6; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (CON)'],[YSTR ' (CON)'],str,cm,useppos(ii,:),20); % D2 CONN
% % % % Summary bar plot
% % %   myfig([YSTR str],useppos(7,:)); B=bar(svcorr,0.6);  hold on; title(TMPLBL2);
% % %   ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %   for pp=1:2:6; B.CData(pp,:)=cmap_pred(1,:); end 
% % %   for pp=2:2:6; B.CData(pp,:)=cmap_pred(2,:); end
% % %   B.FaceColor='flat'; set(gca,'XTick',1.5:2:5.5,opt4{:}); 
% % %   for pp=2:2:4; line([pp+.5 pp+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % % % -- Correlations with BOLD-FC
% % %   D3=FC; YSTR='BOLD-FC'; svcorr=zeros(1,6); useppos=ppostmp3; TMPLBL2=[CMLBL_simple{xx} ' vs ' YSTR];
% % %   X=D1;          Y=D3;          ii=1; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (ALL)'],[YSTR ' (ALL)'],str,cm,useppos(ii,:),20); % D1 all
% % %   X=D2;          Y=D3;          ii=2; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (ALL)'],[YSTR ' (ALL)'],str,cm,useppos(ii,:),20); % D2 all
% % %   X=D1(~masknz); Y=D3(~masknz); ii=3; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (DIS)'],[YSTR ' (DIS)'],str,cm,useppos(ii,:),20); % D1 DIS
% % %   X=D2(~masknz); Y=D3(~masknz); ii=4; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (DIS)'],[YSTR ' (DIS)'],str,cm,useppos(ii,:),20); % D2 DIS
% % %   X=D1(masknz);  Y=D3(masknz);  ii=5; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (CON)'],[YSTR ' (CON)'],str,cm,useppos(ii,:),20); % D1 CONN
% % %   X=D2(masknz);  Y=D3(masknz);  ii=6; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (CON)'],[YSTR ' (CON)'],str,cm,useppos(ii,:),20); % D2 CONN
% % % % Summary bar plot
% % %   myfig([YSTR str],useppos(7,:)); B=bar(svcorr,0.6);  hold on; title(TMPLBL2);
% % %   ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %   for pp=1:2:6; B.CData(pp,:)=cmap_pred(1,:); end 
% % %   for pp=2:2:6; B.CData(pp,:)=cmap_pred(2,:); end
% % %   B.FaceColor='flat'; set(gca,'XTick',1.5:2:5.5,opt4{:}); 
% % %   for pp=2:2:4; line([pp+.5 pp+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % % % Histograms
% % %   myfig(str,ppostmp2(1,:)); histogram(D1,bins); title(Dstr1); font(fntsz,usefont);
% % %   myfig(str,ppostmp2(2,:)); histogram(D2,bins); title(Dstr2); font(fntsz,usefont);
% % % % -- Correlations with each other
% % %   svcorr=zeros(1,3); useppos=ppostmp2; TMPLBL2=['caliber vs myelin (' CMLBL_simple{xx} ')'];
% % %   X=D1;          Y=D2;          [~,svcorr(1)]=ddisp_corr(X,Y,[Dstr1 ' (ALL)'],[Dstr2 ' (ALL)'],str,cm,useppos(3,:),20);
% % %   X=D1(~masknz); Y=D2(~masknz); [~,svcorr(2)]=ddisp_corr(X,Y,[Dstr1 ' (DIS)'],[Dstr2 ' (DIS)'],str,cm,useppos(4,:),20);
% % %   X=D1(masknz);  Y=D2(masknz);  [~,svcorr(3)]=ddisp_corr(X,Y,[Dstr1 ' (CON)'],[Dstr2 ' (CON)'],str,cm,useppos(5,:),20);
% % % % Summary bar plot
% % %   myfig([TMPLBL2 str],useppos(6,:)); B=bar(svcorr,0.6);  hold on; title(TMPLBL2);
% % %   ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %   B.CData=cmaptmp; B.FaceColor='flat'; set(gca,opt4{:}); 
% % % % Connectomes: sequential node ordering
% % %   connplot(D1,[TMPLBL2 str],[],ppostmp(1,:)); title(Dstr1); clim(CLIM1); colormap(cm); set(gca,opt{:})
% % %   connplot(D2,[TMPLBL2 str],[],ppostmp(2,:)); title(Dstr2); clim(CLIM2); colormap(cm); set(gca,opt{:})
% % % % Connectomes: Yeo RSN ordering
% % %   connplot(D1,[TMPLBL2 str],pinfo,ppostmp(3,:)); title(Dstr1); clim(CLIM1); colormap(cm); set(gca,opt2{:})
% % %   connplot(D2,[TMPLBL2 str],pinfo,ppostmp(4,:)); title(Dstr2); clim(CLIM2); colormap(cm); set(gca,opt2{:})
% % % 
% % % 
% % % % ==========   Scatter plot caliber vs myelin: data colored by connectivity
% % % 
% % %   xx=6; 
% % %   tmask=zeros(Nnode); tmask(masknz)=1; ppostmp=[-2558 437 415 360];
% % % 
% % % % ---- zscored data
% % %   D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; % CMAX=prctile(abs(D1(D1~=0)),99); CLIM1=[CMAX*-1 CMAX]; %CMAX=abs(max(D1(:)));
% % %   D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx}; % CMAX=prctile(abs(D2(D2~=0)),99); CLIM2=[CMAX*-1 CMAX]; %CMAX=abs(max(D2(:)));
% % % % -- All edges colored by connectivity
% % %   TMPLBL2=['caliber vs myelin (' CMLBL_simple{xx} ')'];
% % %   X=D1(:); Y=D2(:); tmask_vec=tmask(:);                                     % flatten
% % %   X(maskmd(:))=[]; Y(maskmd(:))=[]; tmask_vec(maskmd(:))=[];                % Remove main diagonal
% % % % Function to compute density and scatter
% % %   plotDensityScatt=@(x,y,baseColor) scatter(x,y,5,densityColor(x,y,baseColor),'filled','MarkerFaceAlpha',0.8);
% % % % Plot
% % %   myfig(TMPLBL2,ppostmp); hold on;
% % %   plotDensityScatt(X(tmask_vec==1), Y(tmask_vec==1), cmap_discon(2,:)); % red = connected
% % %   plotDensityScatt(X(~tmask_vec),Y(~tmask_vec),cmap_discon(1,:)); % blue = disconnected
% % %   xlabel(extractBefore(Dstr1,'-')); ylabel(extractBefore(Dstr2,'-')); title(CMLBL_simple{xx}); font(fntsz,usefont);
% % %   box on;
% % % 
% % % 
% % % % ==========   Scatter plot caliber vs myelin: data colored by ED
% % %   xx=6; 
% % %   ppostmp=[-2558 437 415 360];
% % % % ---- zscored data
% % %   D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; % caliber communication model
% % %   D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx}; % myelin communication model
% % %   D3=EDnrm;
% % % % -- All edges colored by connectivity
% % %   TMPLBL2='caliber vs myelin (data colored by Euclid dist)';
% % %   X=D1(:); Y=D2(:); Z=D3(:);                                                % flatten
% % %   X(maskmd(:))=[]; Y(maskmd(:))=[]; Z(maskmd(:))=[];                        % Remove main diagonal
% % % % Plot
% % %   myfig(TMPLBL2,ppostmp); hold on;
% % %   scatter(X,Y,5,Z,'filled','MarkerFaceAlpha', 0.7);
% % %   colormap(cm3); CB=colorbar; %ylabel(CB,'ED','FontSize',fntsz);
% % %   xlabel(extractBefore(Dstr1,'-')); ylabel(extractBefore(Dstr2,'-')); title(CMLBL_simple{xx}); font(fntsz,usefont);
% % %   box on;
% % %   % ylim([-3 5]); xlim([-2.1 6]); % SPE
% % %   % xlim([-3.9 6]); ylim([-5 5]); % NE
% % %   % ylim([-1.6 3.6]); xlim([-2.2 6]); % SIE
% % %   % xlim([-2.5 4.2]); ylim([-2.5 4.5]); % PT
% % %   % xlim([-3.5 5]); ylim([-3 4]); % DE
% % %   line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','--','LineWidth',3);
% % %   line([0 0],ylim,'Color',[.5 .5 .5],'LineStyle','--','LineWidth',3);
% % % 
% % % 
% % % % ==========   Scatter plot caliber vs myelin: data colored by Yeo RSNs
% % %   xx=6; 
% % %   D=pinfo.masks.rsnetwork; D(isnan(D))=0; D=D+D.'; tmask=D; 
% % %   cmaptmp=[.5 .5 .5; cmap_ntwk]; ppostmp=[-2558 437 415 360];
% % % % ---- zscored data
% % %   D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; % caliber communication model
% % %   D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx}; % myelin communication model
% % % % -- All edges colored by connectivity
% % %   TMPLBL2='caliber vs myelin (data colored by Yeo RSNs)';
% % %   X=D1(:); Y=D2(:); tmask_vec=tmask(:);                                     % flatten
% % %   X(maskmd(:))=[]; Y(maskmd(:))=[]; tmask_vec(maskmd(:))=[];                % Remove main diagonal
% % % % Plot
% % %   myfig(TMPLBL2,ppostmp); hold on;
% % %   for ii=0:Nntwk
% % %     idx=(tmask_vec==ii); scatter(X(idx),Y(idx),10,cmaptmp(ii+1,:),'filled','MarkerFaceAlpha',0.5);
% % %   end
% % %   xlabel(extractBefore(Dstr1,'-')); ylabel(extractBefore(Dstr2,'-')); title(CMLBL_simple{xx}); font(fntsz,usefont);
% % %   box on;
% % %   % ylim([-3 5]); xlim([-2.1 6]); % SPE
% % %   % xlim([-3.9 6]); ylim([-5 5]); % NE
% % %   % ylim([-1.6 3.6]); xlim([-2.2 6]); % SIE
% % %   % xlim([-2.5 4.2]); ylim([-2.5 4.5]); % PT
% % %   xlim([-3.5 5]); ylim([-3 4]); % DE
% % %   line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','--','LineWidth',3);
% % %   line([0 0],ylim,'Color',[.5 .5 .5],'LineStyle','--','LineWidth',3);
% % % 
% % % 
% % % % ==========   BOXPLOT caliber vs myelin: data grouped by Yeo RSNs
% % %   xx=1; 
% % %   D=pinfo.masks.rsnetwork; D(isnan(D))=0; D=D+D.'; tmask=D; 
% % %   cmaptmp=[.5 .5 .5; cmap_ntwk]; tmplbls=['BTW'; LBL_ntwk]; ppostmp=[-2558 482 988 315];
% % %   boxopt={'MarkerStyle','none','BoxWidth',0.35,'BoxFaceAlpha',0.5,'LineWidth',1.5,'BoxEdgeColor','k'};
% % %   scatopt={'filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1};
% % %   textopt={'***','HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',fntsz,'FontWeight','bold'};
% % % % ---- zscored data
% % %   D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx};                                      % caliber communication model
% % %   D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx};                                      % myelin communication model
% % % % -- All edges colored by connectivity
% % %   TMPLBL2='caliber vs myelin (data grouped by Yeo RSNs)';
% % %   X=D1(:); Y=D2(:); tmask_vec=tmask(:);                                     % flatten
% % %   X(maskmd(:))=[]; Y(maskmd(:))=[]; tmask_vec(maskmd(:))=[];                % Remove main diagonal
% % % % Plot
% % %   myfig(TMPLBL2,ppostmp); hold on;
% % %   for ii=0:Nntwk
% % %         idx=(tmask_vec==ii); tx=X(idx); ty=Y(idx); x_tx=ii-0.2; x_ty=ii+0.2;
% % %       % Plot data points first
% % %         jit=0.1;
% % %         scatter(x_tx+(rand(size(tx))-0.5)*2*jit,tx,25,cmap_pred(1,:),scatopt{:});
% % %         scatter(x_ty+(rand(size(ty))-0.5)*2*jit,ty,25,cmap_pred(2,:),scatopt{:});
% % %       % Overlay Box plots
% % %         boxchart(ones(size(tx))*x_tx,tx,'BoxFaceColor',cmap_pred(1,:),boxopt{:}); % Caliber
% % %         boxchart(ones(size(ty))*x_ty,ty,'BoxFaceColor',cmap_pred(2,:),boxopt{:}); % Myelin   
% % %       % Asterisks for t-test
% % %         [~,p]=ttest(tx,ty);
% % %       if p<0.001
% % %         maxy=6; 
% % %         %maxy=max([tx;ty])*1.05;
% % %         plot([x_tx x_ty],[maxy maxy],'k-','LineWidth',1.5);
% % %         text(mean([x_tx x_ty]),maxy,textopt{:});
% % %       end
% % %   end
% % % % Formatting
% % %   xticks(0:Nntwk); xticklabels(tmplbls); %xtickangle(45); xticks(-1.5:3:3*Nntwk);
% % %   ylabel(CMLBL_simple{xx}); font(fntsz,usefont); % legend({extractBefore(Dstr1,'-'),extractBefore(Dstr2,'-')}, 'Location','Best');
% % %   box on; grid on; set(gca,'XGrid','off'); font(fntsz,usefont);
% % %   line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','--','LineWidth',3);
% % % 
% % % 
% % % % ==========   CMs vs all FC
% % % 
% % %   % % % xx=6; 
% % % % % %   ppostmp=[-539 -349 -160 28 218 407 597]; useY=Yz_m([8:-1:3 1]); YLBL=ylbl([8:-1:3 1]);
% % % % % %   opt4={'XTickLabel',{'all' 'dis' 'con'},'XTickLabelRotation',0}; opt={'XTickLabels','','YTickLabels',''};
% % % % % % 
% % % % % % % ---- zscored data
% % % % % %   D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; CMAX=prctile(abs(D1(D1~=0)),99); CLIM1=[CMAX*-1 CMAX];
% % % % % %   D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx}]; CMAX=prctile(abs(D2(D2~=0)),99); CLIM2=[CMAX*-1 CMAX];
% % % % % % % --- Loop over FC
% % % % % %   for yy=1:7
% % % % % %       D3=useY{yy}; Dstr3=YLBL{yy}; CMAX=prctile(abs(D3(D3~=0)),99); CLIM3=[CMAX*-1 CMAX];
% % % % % %       useppos=getPlotPos([-2559 ppostmp(yy) 238 187],10);
% % % % % %     % -- Correlations
% % % % % %       svcorr=zeros(1,6); TMPLBL2=[CMLBL_simple{xx} ' vs ' Dstr3];
% % % % % %       X=D1;          Y=D3;          ii=1; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (ALL)'],[Dstr3 ' (ALL)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end
% % % % % %       X=D2;          Y=D3;          ii=2; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (ALL)'],[Dstr3 ' (ALL)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end; ylabel '';
% % % % % %       X=D1(~masknz); Y=D3(~masknz); ii=3; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (DIS)'],[Dstr3 ' (DIS)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end; ylabel '';
% % % % % %       X=D2(~masknz); Y=D3(~masknz); ii=4; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (DIS)'],[Dstr3 ' (DIS)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end; ylabel '';
% % % % % %       X=D1(masknz);  Y=D3(masknz);  ii=5; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (CON)'],[Dstr3 ' (CON)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end; ylabel '';
% % % % % %       X=D2(masknz);  Y=D3(masknz);  ii=6; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (CON)'],[Dstr3 ' (CON)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end; ylabel '';
% % % % % %     % Summary bar plot
% % % % % %       myfig([Dstr3 str],useppos(7,:)); B=bar(svcorr,0.6);  hold on; % title(TMPLBL2);
% % % % % %       ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % % % % %       for pp=1:2:6; B.CData(pp,:)=cmap_pred(1,:); end 
% % % % % %       for pp=2:2:6; B.CData(pp,:)=cmap_pred(2,:); end
% % % % % %       B.FaceColor='flat'; set(gca,'XTick',1.5:2:5.5,opt4{:}); 
% % % % % %       for pp=2:2:4; line([pp+.5 pp+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % % % % %     % Connectomes: Yeo RSN ordering
% % % % % %       connplot(D1,[TMPLBL2 str],pinfo,useppos(8,:));  title(Dstr1); clim(CLIM1); colormap(cm2); set(gca,opt{:}); font(18,usefont)
% % % % % %       connplot(D2,[TMPLBL2 str],pinfo,useppos(9,:));  title(Dstr2); clim(CLIM2); colormap(cm2); set(gca,opt{:}); font(18,usefont)
% % % % % %       connplot(D3,[TMPLBL2 str],pinfo,useppos(10,:)); title(Dstr3); clim(CLIM3); colormap(cm2); set(gca,opt{:}); font(18,usefont)
% % % % % %   end
% % % 
% % % 
% % % % ==== Box plots of within & between RSNs for CMi caliber vs myelin
% % % % myelin supports greater SPE & NE: between all RSNs, within transmodal
% % % % is myelin-SI also lower in transmodal & between RSNs?
% % % % does PT fit?
% % % % CMY?
% % % % myelin-DE higher in transmodal & between as well?
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % ================ SA-axis vs CM-strength
% % % % scatters & bars
% % % 
% % % % ---- ALL EDGES
% % %   dx=SA_map; dxlbl='S-A map'; savecorr=zeros(2,Ncm); pposvals=[494 111]; str2=' (ALL)';
% % % % -- zscored data
% % %   for ii=1:2
% % %       ppostmp=getPlotPos([-2559 pposvals(ii) 344 303],Ncm+1,0);
% % %       Dstr=[useDstr{ii} str2];
% % %     % Over models
% % %       for xx=1:Ncm
% % %           D=ALLCM{ii}{xx};
% % %           D(isinf(D))=nan; dy=nansum(D,2); CLIM=[min(dy(dy~=0)) max(dy(:))];     
% % %           dylbl=[CMLBL_simple{xx} ' strength'];
% % %         % Scatter plot
% % %           ddisp_corr(dx,dy,dxlbl,dylbl,Dstr,cm,ppostmp(xx,:),fntsz); 
% % %           set(gca,'XTickLabels',''); colorbar off; ylim(CLIM); % ylabel '';
% % %           line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','-','LineWidth',1);
% % %         % Save correlation
% % %           savecorr(ii,xx)=corr(dx(dy~=0),dy(dy~=0),'type','Spearman');
% % %       end
% % %     % Bar plot summary
% % %       D=savecorr(ii,:);
% % %       myfig(Dstr,ppostmp(xx+1,:)); B=bar(CMLBL_simple,D,0.5);
% % %       ylabel('\rho'); font(fntsz,usefont); %title(Dstr);
% % %       grid on; set(gca,'XGrid','off','XMinorGrid','off');
% % %       B.CData=cmap_6cm; B.FaceColor='flat'; % ylim([-.52 .02]);
% % %   end
% % % % -- no xfm data
% % %   for ii=1:2
% % %       ppostmp=getPlotPos([-2559 pposvals(ii) 344 303],Ncm+1,0);
% % %       Dstr=[useDstr{ii} str2];
% % %     % Over models
% % %       for xx=1:Ncm
% % %           D=ALLCMnx{ii}{xx};
% % %           D(isinf(D))=nan; dy=nansum(D,2); CLIM=[min(dy(dy~=0)) max(dy(:))];     
% % %           dylbl=[CMLBL_simple{xx} ' strength'];
% % %         % Scatter plot
% % %           ddisp_corr(dx,dy,dxlbl,dylbl,Dstr,cm,ppostmp(xx,:),fntsz); 
% % %           set(gca,'XTickLabels',''); colorbar off; ylim(CLIM); % ylabel '';
% % %           line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','-','LineWidth',1);
% % %         % Save correlation
% % %           savecorr(ii,xx)=corr(dx(dy~=0),dy(dy~=0),'type','Spearman');
% % %       end
% % %     % Bar plot summary
% % %       D=savecorr(ii,:);
% % %       myfig(Dstr,ppostmp(xx+1,:)); B=bar(CMLBL_simple,D,0.5);
% % %       ylabel('\rho'); font(fntsz,usefont); %title(Dstr);
% % %       grid on; set(gca,'XGrid','off','XMinorGrid','off');
% % %       B.CData=cmap_6cm; B.FaceColor='flat'; % ylim([-.52 .02]);
% % %   end
% % % 
% % % 
% % % % ---- CON EDGES
% % %   dx=SA_map; dxlbl='S-A map'; savecorr=zeros(2,Ncm); pposvals=[494 111]; str2=' (CON)';
% % % % -- zscored data
% % %   for ii=1:2
% % %       ppostmp=getPlotPos([-2559 pposvals(ii) 344 303],Ncm+1,0);
% % %       Dstr=[useDstr{ii} str2];
% % %     % Over models
% % %       for xx=1:Ncm
% % %           D=ALLCM{ii}{xx}; D(~masknz)=0;
% % %           D(isinf(D))=nan; dy=nansum(D,2); CLIM=[min(dy(dy~=0)) max(dy(:))];     
% % %           dylbl=[CMLBL_simple{xx} ' strength'];
% % %         % Scatter plot
% % %           ddisp_corr(dx,dy,dxlbl,dylbl,Dstr,cm,ppostmp(xx,:),fntsz); 
% % %           set(gca,'XTickLabels',''); colorbar off; ylim(CLIM); % ylabel '';
% % %           line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','-','LineWidth',1);
% % %         % Save correlation
% % %           savecorr(ii,xx)=corr(dx(dy~=0),dy(dy~=0),'type','Spearman');
% % %       end
% % %     % Bar plot summary
% % %       D=savecorr(ii,:);
% % %       myfig(Dstr,ppostmp(xx+1,:)); B=bar(CMLBL_simple,D,0.5);
% % %       ylabel('\rho'); font(fntsz,usefont); %title(Dstr);
% % %       grid on; set(gca,'XGrid','off','XMinorGrid','off');
% % %       B.CData=cmap_6cm; B.FaceColor='flat'; % ylim([-.52 .02]);
% % %   end
% % % % -- no xfm data
% % %   for ii=1:2
% % %       ppostmp=getPlotPos([-2559 pposvals(ii) 344 303],Ncm+1,0);
% % %       Dstr=[useDstr{ii} str2];
% % %     % Over models
% % %       for xx=1:Ncm
% % %           D=ALLCMnx{ii}{xx}; D(~masknz)=0;
% % %           D(isinf(D))=nan; dy=nansum(D,2); CLIM=[min(dy(dy~=0)) max(dy(:))];     
% % %           dylbl=[CMLBL_simple{xx} ' strength'];
% % %         % Scatter plot
% % %           ddisp_corr(dx,dy,dxlbl,dylbl,Dstr,cm,ppostmp(xx,:),fntsz); 
% % %           set(gca,'XTickLabels',''); colorbar off; ylim(CLIM); % ylabel '';
% % %           line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','-','LineWidth',1);
% % %         % Save correlation
% % %           savecorr(ii,xx)=corr(dx(dy~=0),dy(dy~=0),'type','Spearman');
% % %       end
% % %     % Bar plot summary
% % %       D=savecorr(ii,:);
% % %       myfig(Dstr,ppostmp(xx+1,:)); B=bar(CMLBL_simple,D,0.5);
% % %       ylabel('\rho'); font(fntsz,usefont); %title(Dstr);
% % %       grid on; set(gca,'XGrid','off','XMinorGrid','off');
% % %       B.CData=cmap_6cm; B.FaceColor='flat'; % ylim([-.52 .02]);
% % %   end
% % % 
% % % % ---- DIS EDGES
% % %   dx=SA_map; dxlbl='S-A map'; savecorr=zeros(2,Ncm); pposvals=[494 111]; str2=' (DIS)';
% % % % -- zscored data
% % %   for ii=1:2
% % %       ppostmp=getPlotPos([-2559 pposvals(ii) 344 303],Ncm+1,0);
% % %       Dstr=[useDstr{ii} str2];
% % %     % Over models
% % %       for xx=1:Ncm
% % %           D=ALLCM{ii}{xx}; D(masknz)=0;
% % %           D(isinf(D))=nan; dy=nansum(D,2); CLIM=[min(dy(dy~=0)) max(dy(:))];     
% % %           dylbl=[CMLBL_simple{xx} ' strength'];
% % %         % Scatter plot
% % %           ddisp_corr(dx,dy,dxlbl,dylbl,Dstr,cm,ppostmp(xx,:),fntsz); 
% % %           set(gca,'XTickLabels',''); colorbar off; ylim(CLIM); % ylabel '';
% % %           line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','-','LineWidth',1);
% % %         % Save correlation
% % %           savecorr(ii,xx)=corr(dx(dy~=0),dy(dy~=0),'type','Spearman');
% % %       end
% % %     % Bar plot summary
% % %       D=savecorr(ii,:);
% % %       myfig(Dstr,ppostmp(xx+1,:)); B=bar(CMLBL_simple,D,0.5);
% % %       ylabel('\rho'); font(fntsz,usefont); %title(Dstr);
% % %       grid on; set(gca,'XGrid','off','XMinorGrid','off');
% % %       B.CData=cmap_6cm; B.FaceColor='flat'; % ylim([-.52 .02]);
% % %   end
% % % % -- no xfm data
% % %   for ii=1:2
% % %       ppostmp=getPlotPos([-2559 pposvals(ii) 344 303],Ncm+1,0);
% % %       Dstr=[useDstr{ii} str2];
% % %     % Over models
% % %       for xx=1:Ncm
% % %           D=ALLCMnx{ii}{xx}; D(masknz)=0;
% % %           D(isinf(D))=nan; dy=nansum(D,2); CLIM=[min(dy(dy~=0)) max(dy(:))];     
% % %           dylbl=[CMLBL_simple{xx} ' strength'];
% % %         % Scatter plot
% % %           ddisp_corr(dx,dy,dxlbl,dylbl,Dstr,cm,ppostmp(xx,:),fntsz); 
% % %           set(gca,'XTickLabels',''); colorbar off; ylim(CLIM); % ylabel '';
% % %           line(xlim,[0 0],'Color',[.5 .5 .5],'LineStyle','-','LineWidth',1);
% % %         % Save correlation
% % %           savecorr(ii,xx)=corr(dx(dy~=0),dy(dy~=0),'type','Spearman');
% % %       end
% % %     % Bar plot summary
% % %       D=savecorr(ii,:);
% % %       myfig(Dstr,ppostmp(xx+1,:)); B=bar(CMLBL_simple,D,0.5);
% % %       ylabel('\rho'); font(fntsz,usefont); %title(Dstr);
% % %       grid on; set(gca,'XGrid','off','XMinorGrid','off');
% % %       B.CData=cmap_6cm; B.FaceColor='flat'; % ylim([-.52 .02]);
% % %   end
% % % 
% % % 
% % % % ================ SA-axis vs CM-strength
% % % % ----- Matrix of ALL CON & DIS together
% % %   useD=ALLCM;     str2=' (z-scored)';                           % zscored data
% % %   % useD=ALLCMnx;   str2=' (no xfm)';                             % no xfm data
% % %   dx=SA_map; dxlbl='S-A map'; savecorr=zeros(2*3,Ncm); 
% % %   TMPLBLS={'ALL' 'DIS' 'CON'};
% % % % -- Compute correlations
% % %   for ii=1:2
% % %     % ALL
% % %       for xx=1:Ncm
% % %           D=useD{ii}{xx}; D(isinf(D))=nan; dy=nansum(D,2); 
% % %           savecorr(ii,xx)=corr(dx(dy~=0),dy(dy~=0),'type','Spearman');
% % %       end
% % %     % DIS
% % %       for xx=1:Ncm
% % %           D=useD{ii}{xx}; D(masknz)=0; D(isinf(D))=nan; dy=nansum(D,2); 
% % %           savecorr(ii+2,xx)=corr(dx(dy~=0),dy(dy~=0),'type','Spearman');
% % %       end
% % %     % CON
% % %       for xx=1:Ncm
% % %           D=useD{ii}{xx}; D(~masknz)=0; D(isinf(D))=nan; dy=nansum(D,2); 
% % %           savecorr(ii+4,xx)=corr(dx(dy~=0),dy(dy~=0),'type','Spearman');
% % %       end
% % %   end
% % % % -- Plot
% % %   ppostmp=getPlotPos([-2559 389 573 408],2,0);
% % %   D=savecorr; D=D([1 3 5 2 4 6],:); D=D';
% % %   TMPSTR1='SA-axis vs CM-strength'; TMPSTR2=[TMPSTR1 str2];
% % %   CLIM=[-.6 .6]; %CMAX=ceil(max(abs(D(D~=0)))*10)/10; CLIM=[CMAX*-1 CMAX];
% % %   myfig([TMPSTR2 str],ppostmp(1,:)); imagesc(D); CB=colorbar; colormap(cm2); 
% % %   clim(CLIM); font(fntsz,usefont); % title(TMPSTR2); 
% % %   set(gca,'YTickLabels',CMLBL_simple,'XTick',[1:6],'XTickLabels',[TMPLBLS TMPLBLS],'XTickLabelRotation',45);
% % %   ylabel(CB,'\rho','FontSize',fntsz,'FontWeight','bold','FontName',usefont,'Rotation',0);
% % % % - Add labels to top of x-axis
% % % % Duplicate axes
% % %   ax1=gca; ax2=axes('Position',ax1.Position,'XAxisLocation','top','YAxisLocation','right','Color','none','XColor','k', 'YColor','none');
% % % % Match limits
% % %   ax2.XLim=ax1.XLim;  ax2.YLim=ax1.YLim;
% % % % Add labels to top axis
% % %   set(ax2,'XTick',[2 5],'XTickLabel',{'caliber' 'myelin'},'TickLength',[0 0]);
% % %   set(ax2, 'Box','off','FontWeight','bold');
% % %   font(fntsz,usefont);
% % % 
% % % % -------------------------------------------------------------------------
% % % 
% % % 
% % % 
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % %% =============      Modeling 1a: FCx ~ CMi_Wj + ED        ============ %
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % FCx modeled by a single weighted communication model.
% % % % Caliber vs myelin
% % % 
% % % Nw8s=2; Ncoeff=2; 
% % % PLTLBL_hi=cell(1,Nw8s);
% % % 
% % % 
% % % % -------------------------------------------------------------------------
% % % %% ==========================      Global      ============================
% % % % -------------------------------------------------------------------------
% % % 
% % %   CC=2;                    % Switch for communication models
% % %   SNP=0;                   % switch for ALL (0) vs CON (1) vs DIS (2) edges
% % % 
% % % % Euclidean distance
% % %   Xd=EDnrm; DISTLBL='ED'; Xdz=nzzscore(Xd,0);
% % % 
% % % % GLOBAL model variables
% % %   R2_g=zeros(Nw8s,Nfc);                RMSE_g=zeros(Nw8s,Nfc);   
% % %   PVAL_g=zeros(Ncoeff,Nfc,Nw8s);      
% % %   DA_g_value=zeros(Ncoeff,Nfc,Nw8s);   DA_g_prcnt=zeros(Ncoeff,Nfc,Nw8s);  
% % %   DA_g_ranks=zeros(Ncoeff,Nfc,Nw8s);   B_g=zeros(Ncoeff,Nfc,Nw8s); 
% % %   y_pred_g=cell(Nw8s,Nfc);             corr_pred_g=zeros(Nw8s,Nfc);
% % % 
% % % % NETWORK model variables
% % %   R2_ntwk=zeros(Nntwk,Nntwk,Nfc,Nw8s);          DA_ranks_ntwk=zeros(Nntwk,Nntwk,Ncoeff,Nfc,Nw8s);
% % %   DFE_ntwk=zeros(Nntwk,Nntwk,Nfc,Nw8s);         DA_value_ntwk=zeros(Nntwk,Nntwk,Ncoeff,Nfc,Nw8s);
% % %   RMSE_ntwk=zeros(Nntwk,Nntwk,Nfc,Nw8s);             DA_prcnt_ntwk=zeros(Nntwk,Nntwk,Ncoeff,Nfc,Nw8s);         
% % %   B_ntwk=zeros(Nntwk,Nntwk,Ncoeff,Nfc,Nw8s);    PVAL_ntwk=zeros(Nntwk,Nntwk,Ncoeff,Nfc,Nw8s);  
% % %   y_pred_ntwk=cell(Nntwk,Nntwk,Nfc,Nw8s);       corr_pred_ntwk=zeros(Nntwk,Nntwk,Nfc,Nw8s);
% % %   corr_pred_ntwk_pool=zeros(Nw8s,Nfc);
% % % 
% % % % NODE model variables
% % %   B_node=zeros(Ncoeff,Nnode,Nfc,Nw8s);          R2_node=zeros(Nnode,Nfc,Nw8s); 
% % %   PVAL_node=zeros(Ncoeff,Nnode,Nfc,Nw8s);       DFE_node=zeros(Nnode,Nfc,Nw8s); 
% % %   DA_node_value=zeros(Ncoeff,Nnode,Nfc,Nw8s);   DA_node_prcnt=zeros(Ncoeff,Nnode,Nfc,Nw8s);
% % %   DA_node_ranks=zeros(Ncoeff,Nnode,Nfc,Nw8s);   
% % %   y_pred_node=cell(Nnode,Nfc,Nw8s);             corr_pred_node=zeros(Nnode,Nfc,Nw8s);
% % %   PLTLBL_hi=cell(1,Nw8s);                       corr_pred_node_pool=zeros(Nw8s,Nfc);
% % % 
% % % 
% % % 
% % % % ===== Loop over caliber & myelin
% % %   for jj=1:2
% % % 
% % %     % Communication model for caliber or myelin
% % %       Xz=ALLCM{jj}{CC}; Xnx=ALLCMnx{jj}{CC}; Xstr=[CMLBLnx{jj}{CC}]; Xstr=strrep(Xstr,'-','');
% % % 
% % %     % Combine all X predictors
% % %       useD={Xz  Xdz};     useDnx={Xnx  Xd};    DLBL={Xstr  DISTLBL};
% % %       useN=length(useD); ixd=useN; nd=ixd+1; iy=nd; PLTLBL=DLBL;
% % % 
% % % 
% % %     % Filter for ALL vs CON vs DIS edges
% % %       switch SNP 
% % %           case 0
% % %               Dcdastr=' (ALL-NP)';
% % %           case 1
% % %               for kk=1:ixd;  useD{kk}=useD{kk}(masknz); end
% % %               Dcdastr=' (CON-NP)';
% % %           case 2
% % %               for kk=1:ixd;  useD{kk}=useD{kk}(~masknz); end
% % %               Dcdastr=' (DIS-NP)';
% % %       end
% % % 
% % % 
% % % 
% % %     % ---- Combine & compute metadata
% % %     % Combine X & Y
% % %       ALLD=cell(1,nd);   for ii=1:useN; ALLD{ii}=useD{ii};     end; ALLD{iy}=Yz_m; 
% % %       ALLDnx=cell(1,nd); for ii=1:useN; ALLDnx{ii}=useDnx{ii}; end; ALLDnx{iy}=Y_m;
% % % 
% % %     % ---- MODEL SPECS
% % %     % Model specs 1: FC ~ 1 + CMi_Wj + ED
% % %       mdlstrgen='y~1+x1+x2';
% % %       mdlstrlonggen=mdlstrgen;
% % %       pltlbl_mdl=['FC ~ 1 + ' PLTLBL{1} ' + ' PLTLBL{2}];  
% % %       D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{1};
% % %       Dt=rmzeros(D);  X=Dt(:,1:ixd); y=Dt(:,iy);                              % Remove 0-valued edges
% % %     % Modeling
% % %       mdlspec=[ylbl{1} '~1+' PLTLBL{1} '+' PLTLBL{2}];
% % %       M=fitlm(X,y,mdlspec,'VarNames',[PLTLBL ylbl{1}]);
% % %       PLTLBL_hi{jj}=M.Coefficients.Properties.RowNames(2:end)';
% % % 
% % % 
% % %     % ---- Global modeling
% % %       for ii = 1 : Nfc
% % %         % Extract all X data with the FC data for this round
% % %           D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
% % %           % Dnx=cell(1,iy); Dnx(1:ixd)=ALLDnx(1:ixd); Dnx{iy}=ALLDnx{iy}{ii};
% % %         % Check for presence of subcortical data in FC
% % %           if ~all(cellfun(@(c) length(c)==Nnode,D))
% % %               D(1:ixd)=cellfun(@(c) c(Nstx+1:end,Nstx+1:end), D(1:ixd),'UniformOutput',0);
% % %           end
% % %         % Remove 0-valued edges
% % %           Dt=rmzeros(D); X=Dt(:,1:ixd); y=Dt(:,iy); % Dtnx=rmzeros(Dnx);
% % %         % Specify model
% % %           mdlspec=[ylbl{ii} '~1+' PLTLBL{1} '+' PLTLBL{2}]; 
% % %           mdlstrfull=mdlspec; 
% % %         % Modeling
% % %           M=fitlm(X,y,mdlspec,'VarNames',[PLTLBL ylbl{ii}]); 
% % %         % Dominance analysis
% % %           [tmp,~]=dominance_interaction(X,y,mdlstrgen); % DA_g_hi.Variable_name=PLTLBL_hi';
% % %           DA_g_prcnt(:,ii,jj)=tmp.Dominance_percentage;                              
% % %           DA_g_value(:,ii,jj)=tmp.Dominance_value;
% % %           DA_g_ranks(:,ii,jj)=tmp.Rank;
% % %         % Extract model coeffs
% % %           B_g(:,ii,jj)=M.Coefficients.Estimate(2:end);  R2_g(jj,ii)=M.Rsquared.Adjusted; 
% % %           PVAL_g(:,ii,jj)=M.Coefficients.pValue(2:end); RMSE_g(jj,ii)=M.RMSE;
% % %         % Store empirical & predicted y values
% % %           y_pred_g{jj,ii}=[y M.Fitted]; corr_pred_g(jj,ii)=corr(y, M.Fitted);
% % %       end
% % % 
% % %   end
% % % 
% % % 
% % % % % % % ---- Plot model outputs
% % % % % %   str2=['Global model' Dcdastr]; 
% % % % % % % -- Combine R2 & dominance
% % % % % %   YLIM=ceil(max(R2_g(:))*10)/10; ppostmp=[-2559 405 429 392; -2560 -67 429 392];
% % % % % %   for jj=1:2
% % % % % %         D1=R2_g(jj,:); D2=DA_g_prcnt(:,:,jj) .* .01; DPLT=D1.*D2;
% % % % % %         myfig([str2 str],ppostmp(jj,:)); B=bar(ylbl,DPLT',0.5,'stacked'); font(fntsz,usefont);
% % % % % %         ylabel('adjusted R^2'); title(strjoin(PLTLBL_hi{jj}, ' + ')); grid on; set(gca,'XGrid','off');
% % % % % %         B(1).FaceColor=cmap_pred(jj,:); B(2).FaceColor=cmap_pred(4,:); ylim([0 YLIM]);
% % % % % %   end
% % % % % % 
% % % % % % % -- Betas (bar plot + LINES)
% % % % % %   ppostmp=[-2129 405 1209 392; -2130 -67 1209 392];
% % % % % %   D=B_g; D(PVAL_g > 0.05)=0; YLIM=[floor(min(D(:))*10)/10  ceil(max(D(:))*10)/10];
% % % % % %   for jj=1:2
% % % % % %       DPLT=D(:,:,jj);
% % % % % %       myfig([str2 str],ppostmp(jj,:)); B=bar(DPLT',0.8); hold on; font(fntsz,usefont);
% % % % % %       ylabel('\beta p<0.05'); grid on;  title(strjoin(PLTLBL_hi{jj}, ' + '));
% % % % % %       B(1).FaceColor=cmap_pred(jj,:); B(2).FaceColor=cmap_pred(4,:); ylim(YLIM);
% % % % % %       for ii=1:Nfc-1; line([ii+.5 ii+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % % % % %       set(gca,'XTickLabels',ylbl);
% % % % % %   end
% % % % % % 
% % % % % % % -- Empirical vs Predicted
% % % % % % % Scatter plots
% % % % % %   % % % ppostmp=getPlotPos([-2559 359 360 330],Nfc,1);
% % % % % %   % % % for jj=1:Nw8s
% % % % % %   % % %     for ii=1%:Nfc
% % % % % %   % % %         Dt=y_pred_g{jj,ii}; Dx=Dt(:,1); Dy=Dt(:,2); 
% % % % % %   % % %         myfig([PLTLBL_hi{jj}{1} '; ' str2],ppostmp(ii,:)); heatscatter(Dx,Dy); colormap(cm);
% % % % % %   % % %         xlabel('empirical'); ylabel('predicted'); font(fntsz,usefont); 
% % % % % %   % % %         title([sprintf('r=%.2f (',corr_pred_g(jj,ii)) ylbl{ii} ')']);
% % % % % %   % % %     end
% % % % % %   % % % end
% % % % % % % Bar plot
% % % % % %   ppostmp=[-919 405 467 392; -919 -67 467 392];
% % % % % %   for jj=1:Nw8s
% % % % % %       myfig([PLTLBL_hi{jj}{1} '; ' str2],ppostmp(jj,:)); B=bar(ylbl,corr_pred_g(jj,:),0.5); 
% % % % % %       ylabel('Pearson''s r'); title('empirical vs predicted'); font(fntsz,usefont); 
% % % % % %       B.CData=cmap_pred(jj,:); B.FaceColor='flat'; ylim([0 1]); grid on; set(gca,'XGrid','off');
% % % % % %   end
% % % 
% % % 
% % % 
% % % % -------------------------------------------------------------------------
% % % %% =======================    Pairwise networks    ========================
% % % % -------------------------------------------------------------------------
% % % 
% % % % ===== Loop over caliber & myelin
% % %   for jj=1:2
% % % 
% % %     % Communication model for caliber or myelin
% % %       Xz=ALLCM{jj}{CC}; Xstr=CMLBLnx{jj}{CC}; % Xnx=ALLCMnx{jj}{CC};
% % % 
% % %     % Combine all X predictors
% % %       useD={Xz  Xdz};     DLBL={Xstr  DISTLBL}; % useDnx={Xnx  Xd};
% % %       useN=length(useD); ixd=useN; nd=ixd+1; iy=nd; PLTLBL=DLBL;
% % % 
% % % 
% % %     % Filter for ALL vs CON vs DIS edges
% % %       switch SNP 
% % %           case 0
% % %               Dcdastr=' (ALL-NP)';
% % %           case 1
% % %               for kk=1:ixd;  useD{kk}=useD{kk}(masknz); end
% % %               Dcdastr=' (CON-NP)';
% % %           case 2
% % %               for kk=1:ixd;  useD{kk}=useD{kk}(~masknz); end
% % %               Dcdastr=' (DIS-NP)';
% % %       end
% % % 
% % % 
% % %     % ---- Combine & compute metadata        
% % %     % Combine X & Y
% % %       ALLD=cell(1,nd);   for ii=1:useN; ALLD{ii}=useD{ii};     end; ALLD{iy}=Yz_m; 
% % %       % ALLDnx=cell(1,nd); for ii=1:useN; ALLDnx{ii}=useDnx{ii}; end; ALLDnx{iy}=Y_m;
% % % 
% % %     % ---- MODEL SPECS
% % %     % Model specs 1: FC ~ 1 + CMi_Wj + ED
% % %       mdlstrgen='y~1+x1+x2';
% % %       mdlstrlonggen=mdlstrgen;
% % %       pltlbl_mdl=['FC ~ 1 + ' PLTLBL{1} ' + ' PLTLBL{2}];  
% % %       D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{1};
% % %       Dt=rmzeros(D);  X=Dt(:,1:ixd); y=Dt(:,iy);                              % Remove 0-valued edges
% % %     % Modeling
% % %       mdlspec=[ylbl{1} '~1+' PLTLBL{1} '+' PLTLBL{2}];
% % %       M=fitlm(X,y,mdlspec,'VarNames',[PLTLBL ylbl{1}]);
% % %       PLTLBL_hi{jj}=M.Coefficients.Properties.RowNames(2:end)';
% % %       % % % mdlstrlong=strjoin(PLTLBL_hi,' + ');
% % %       % % % mdlstrlong_compact=strjoin(PLTLBL_hi,'+');
% % %       % % % Ncoeff=size(M.Coefficients,1)-1;
% % %       % % % i_intterms=cellstrfind(PLTLBL_hi,':');
% % % 
% % %     % ---- Modeling
% % % 
% % %       for ii = 1 : Nfc
% % %         % Extract all X data with the FC data for this round
% % %           D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
% % %         % Cut the subcortical data out of the SC if not present in FC
% % %           it=0; usepinfo=pinfo;
% % %           if ~all(cellfun(@(c) length(c)==Nnode,D))
% % %               it=1; usepinfo=pinfo2;
% % %               D(1:ixd)=cellfun(@(c) c(Nstx+1:end,Nstx+1:end), D(1:ixd),'UniformOutput',0);
% % %           end
% % %           useNntwk=length(usepinfo.clabels);
% % %         % Specify model
% % %           mdlspec=[ylbl{ii} '~1+' PLTLBL{1} '+' PLTLBL{2}]; 
% % %           disp(['Running model for y-data: ' ylbl{ii}])
% % % 
% % %           for iio = 1 : useNntwk
% % %               nio=usepinfo.cis{iio};
% % %               for iii = iio : useNntwk
% % %                   nii=usepinfo.cis{iii};
% % %                   ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(nio,nii); end 
% % %                   Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
% % %                   try  
% % %                       M=fitlm(X,y,mdlspec,'VarNames',[PLTLBL ylbl{ii}]); 
% % %                       if M.DFE <= 9; error(['Insufficient DFE: ' usepinfo.clabels{iio} '-' usepinfo.clabels{iii}]); end
% % %                       [tmp,~]=dominance_interaction(X,y,mdlstrgen); 
% % %                       DA_value_ntwk(iio+it,iii+it,:,ii,jj)=tmp.Dominance_value;      R2_ntwk(iio+it,iii+it,ii,jj)=M.Rsquared.Adjusted;
% % %                       DA_prcnt_ntwk(iio+it,iii+it,:,ii,jj)=tmp.Dominance_percentage; PVAL_ntwk(iio+it,iii+it,:,ii,jj)=M.Coefficients.pValue(2:end);
% % %                       DA_ranks_ntwk(iio+it,iii+it,:,ii,jj)=tmp.Rank;                 B_ntwk(iio+it,iii+it,:,ii,jj)=M.Coefficients.Estimate(2:end); 
% % %                       DFE_ntwk(iio+it,iii+it,ii,jj)=M.DFE;                           RMSE_ntwk(iio+it,iii+it,ii,jj)=M.RMSE;
% % %                     % Store empirical & predicted y values
% % %                       y_pred_ntwk{iio+it,iii+it,ii,jj}=[y M.Fitted];                   corr_pred_ntwk(iio+it,iii+it,ii,jj)=corr(y,M.Fitted);
% % %                   catch
% % %                       disp(['FAILED BLOCK: Network-pair ' usepinfo.clabels{iio} '-' usepinfo.clabels{iii}])
% % %                       DA_value_ntwk(iio+it,iii+it,:,ii,jj)=zeros(Ncoeff,1);   R2_ntwk(iio+it,iii+it,ii,jj)=0;               
% % %                       DA_ranks_ntwk(iio+it,iii+it,:,ii,jj)=zeros(Ncoeff,1);   B_ntwk(iio+it,iii+it,:,ii,jj)=zeros(Ncoeff,1);              
% % %                       DA_prcnt_ntwk(iio+it,iii+it,:,ii,jj)=zeros(Ncoeff,1);   RMSE_ntwk(iio+it,iii+it,ii,jj)=inf; 
% % %                       DFE_ntwk(iio+it,iii+it,ii,jj)=0;                        PVAL_ntwk(iio+it,iii+it,:,ii,jj)=zeros(Ncoeff,1);
% % % 
% % %                   end
% % %               end
% % %           end
% % %       end
% % % 
% % %     % Pool Pred & Emp y for all nodes
% % %       for ii=1:Nfc
% % %           Dt=y_pred_ntwk(:,:,ii,jj); Dt=Dt(:); Dt(cellfun(@(c) isempty(c),Dt))=[];
% % %           Dx=cell2mat(cellfun(@(c) c(:,1),Dt,'UniformOutput',0));
% % %           Dy=cell2mat(cellfun(@(c) c(:,2),Dt,'UniformOutput',0));
% % %           corr_pred_ntwk_pool(jj,ii)=corr(Dx,Dy,'type','Pearson');
% % %       end
% % %   end
% % % 
% % % 
% % % % ---- Plotting Network-Level (pairwise)
% % % % % % % Pairwise Networks
% % % % % %   str2=['Network-level' Dcdastr];
% % % % % % 
% % % % % % % --- rsNetworks on surface
% % % % % % % Surface plot
% % % % % %   D=[1:Nntwk]'; useN=1; TMPLBL='RSNs'; TMPSTR=['Yeo 7-rsNetworks: ' str]; usecmap=[ .75 .75 .75 ; cmap_ntwk];
% % % % % %   D2=nan(Nnode,useN); for ii=1:Nntwk; for jj=1:useN; D2(pinfo.cis{ii},jj)=D(ii,jj); end; end
% % % % % %   Hcx=plot_conn_surf(D2,pinfo,'cortex',TMPLBL, TMPSTR); setsurf(Hcx,[0 Nntwk],usecmap);
% % % % % % % Legend
% % % % % %   myfig; for ii=1:Nntwk; plot([ii ii],1:2,'-','LineWidth',10,'Color',cmap_ntwk(ii,:)); hold on; end
% % % % % %   legend(LBL_ntwk); font(fntsz,usefont);
% % % % % % 
% % % % % % 
% % % % % % % ----------- Matrix Plots
% % % % % %   pposvals=[498  119];
% % % % % %   axopt={'XTick',1:Nntwk,'XTickLabel',LBL_ntwk,'XTickLabelRotation',40, ...
% % % % % %          'YTick',1:Nntwk,'YTickLabel',LBL_ntwk,'YTickLabelRotation',40,'FontSize',fntsz};
% % % % % %   meshopt={'EdgeColor','k','LineWidth',1.5};
% % % % % % % --- R2
% % % % % %   D=R2_ntwk; D(D<0)=0; TMPSTR='R^2: '; usecmap=cm3; CLIM=[0 ceil(max(D(:))*10)/10];
% % % % % %   for jj=1:Nw8s
% % % % % %       TMPSTR2=[TMPSTR PLTLBL_hi{jj}{1}]; ppostmp2=getPlotPos([-2559 pposvals(jj) 315 299],Nfc,0);
% % % % % %       for ii=1:Nfc
% % % % % %         % Matrix with blank upper tri
% % % % % %           D2=D(:,:,ii,jj)'; tmask=tril(true(size(D2))); myfig([TMPSTR2 ', ' ylbl{ii}],ppostmp2(ii,:));
% % % % % %           imagesc(D2,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; clim(CLIM); font(fntsz,usefont);
% % % % % %           set(gca,axopt{:}); font(fntsz,usefont); hold on; % title(ylbl{ii}); 
% % % % % %           for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % % % % %           set(gca,'XTickLabel','','YTickLabel',''); title ''; colorbar off; axis off;
% % % % % %       end
% % % % % %   end
% % % % % % % --- Beta (p<0.05)
% % % % % %   D=B_ntwk; D(PVAL_ntwk>0.05)=0; TMPSTR='beta: '; usecmap=cm2;  
% % % % % %   for jj=1:Nw8s
% % % % % %       Dt=squeeze(D(:,:,1,:,jj)); CMAX=round(prctile(abs(Dt(:)),99),2); CLIM=[CMAX*-1  CMAX];
% % % % % %       TMPSTR2=[TMPSTR PLTLBL_hi{jj}{1}]; ppostmp2=getPlotPos([-2559 pposvals(jj) 315 299],Nfc,0);
% % % % % %       for ii=1:Nfc
% % % % % %           D2=Dt(:,:,ii)'; tmask=tril(true(size(D2))); myfig([TMPSTR2 ', '  ylbl{ii}],ppostmp2(ii,:));
% % % % % %           imagesc(D2,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; clim(CLIM); font(fntsz,usefont);
% % % % % %           set(gca,axopt{:}); font(fntsz,usefont); hold on; % title(ylbl{ii}); 
% % % % % %           for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % % % % %           set(gca,'XTickLabel','','YTickLabel',''); title ''; axis off; % colorbar off;
% % % % % %       end
% % % % % %   end
% % % % % % % Beta (p<0.05) MEDIAN across FC 
% % % % % %   D=B_ntwk; D(PVAL_ntwk>0.05)=0; TMPSTR='Median \beta_p_<_0_._0_5: '; usecmap=cm2;
% % % % % %   ppostmp=getPlotPos([-2559 446 375 351],2);
% % % % % %   for jj=1:Nw8s 
% % % % % %       Dt=squeeze(D(:,:,1,:,jj)); TMPSTR2=[TMPSTR PLTLBL_hi{jj}{1}]; 
% % % % % %       Dt(Dt==0)=nan; Dt=nanmedian(Dt,3); Dt(isnan(Dt))=0; Dt=Dt'; %CLIM=[-1 1];
% % % % % %       CMAX=round(prctile(abs(Dt(:)),99),2); CLIM=[CMAX*-1  CMAX];
% % % % % %       tmask=tril(true(size(Dt))); myfig(TMPSTR2,ppostmp(jj,:)); 
% % % % % %       imagesc(Dt,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; font(fntsz,usefont); 
% % % % % %       set(gca,axopt{:}); font(fntsz,usefont); hold on; clim(CLIM);
% % % % % %       for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % % % % %       set(gca,'XTickLabel','','YTickLabel',''); title ''; axis off; %colorbar off;
% % % % % %   end
% % % % % % 
% % % % % % % --- Dominance rankings
% % % % % %   D=DA_ranks_ntwk; D(D==0)=nan; D=(D.*-1)+Ncoeff+1; D(isnan(D))=0; 
% % % % % %   TMPSTR='Dominance Rank: '; usecmap=cm4;  CLIM=[0 Ncoeff];
% % % % % %   for jj=1:Nw8s
% % % % % %       Dt=squeeze(D(:,:,1,:,jj)); TMPSTR2=[TMPSTR PLTLBL_hi{jj}{1}]; ppostmp2=getPlotPos([-2559 pposvals(jj) 315 299],Nfc,0); 
% % % % % %       for ii=1:Nfc
% % % % % %           D2=Dt(:,:,ii)'; tmask=tril(true(size(D2))); myfig([TMPSTR2 ', ' ylbl{ii}],ppostmp2(ii,:));
% % % % % %           imagesc(D2,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; clim(CLIM); font(fntsz,usefont);
% % % % % %           set(gca,axopt{:}); font(fntsz,usefont); hold on; % title(ylbl{ii}); 
% % % % % %           for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % % % % %           set(gca,'XTickLabel','','YTickLabel',''); title ''; axis off; %colorbar off;
% % % % % %       end
% % % % % %   end
% % % % % % % Dominance rank MEDIAN across FC
% % % % % %   D=DA_ranks_ntwk; D(D==0)=nan; D=(D.*-1)+Ncoeff+1; usecmap=cm4; CLIM=[0 Ncoeff];
% % % % % %   ppostmp=getPlotPos([-2559 446 375 351],2); TMPSTR='Median dominance rank: ';
% % % % % %   for jj=1:Nw8s 
% % % % % %       D2=squeeze(D(:,:,1,:,jj));  TMPSTR2=[TMPSTR PLTLBL_hi{jj}{1}];
% % % % % %       DPLT=nanmedian(D2,3); DPLT(isnan(DPLT))=0; DPLT=DPLT'; 
% % % % % %       TMPSTR2=[TMPSTR PLTLBL_hi{jj}{1}];  
% % % % % %       tmask=tril(true(size(DPLT))); myfig(TMPSTR2,ppostmp(jj,:));
% % % % % %       imagesc(DPLT,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; clim(CLIM); font(fntsz,usefont);
% % % % % %       set(gca,axopt{:}); font(fntsz,usefont); hold on;
% % % % % %       for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % % % % %       set(gca,'XTickLabel','','YTickLabel',''); title ''; colorbar off; axis off;
% % % % % %   end
% % % % % % 
% % % % % % 
% % % % % % % ----- Empirical vs Predicted:
% % % % % % % Scatter plot: all edges pooled
% % % % % %   % % % pposvals=[498  119];
% % % % % %   % % % for jj=1:Nw8s
% % % % % %   % % %     ppostmp2=getPlotPos([-2559 pposvals(jj) 315 299],Nfc,0);
% % % % % %   % % %     for ii=1:Nfc
% % % % % %   % % %         Dt=y_pred_ntwk(:,:,ii,jj); Dt=Dt(:); Dt(cellfun(@(c) isempty(c),Dt))=[]; 
% % % % % %   % % %         Dx=cell2mat(cellfun(@(c) c(:,1),Dt,'UniformOutput',0));
% % % % % %   % % %         Dy=cell2mat(cellfun(@(c) c(:,2),Dt,'UniformOutput',0));
% % % % % %   % % %         myfig([PLTLBL_hi{jj}{1} str2],ppostmp2(ii,:)); heatscatter(Dx,Dy); colormap(gray);
% % % % % %   % % %         xlabel('empirical'); ylabel('predicted'); font(fntsz,usefont); 
% % % % % %   % % %         title([sprintf('r=%.2f (',corr_pred_ntwk_pool(jj,ii)) ylbl{ii} ')']);
% % % % % %   % % %         colorbar off;
% % % % % %   % % %     end
% % % % % %   % % % end
% % % % % % % Bar plot: edges pooled
% % % % % %   ppostmp=getPlotPos([-2559 467 467 330],2,0);
% % % % % %   for jj=1:Nw8s
% % % % % %       myfig([PLTLBL_hi{jj}{1} ', ' str2],ppostmp(jj,:)); B=bar(ylbl,corr_pred_ntwk_pool(jj,:),0.5); 
% % % % % %       ylabel('Pearson''s r'); title('empirical vs predicted'); font(fntsz,usefont); grid on;
% % % % % %       B.CData=cmap_pred(jj,:); B.FaceColor='flat'; ylim([0 1]);
% % % % % %   end
% % % 
% % % 
% % % 
% % % % -------------------------------------------------------------------------
% % % %% =======================        Node-level       ========================
% % % % -------------------------------------------------------------------------
% % % 
% % % 
% % % % ===== Loop over caliber & myelin
% % %   for jj=1:2
% % % 
% % %     % Communication model for caliber or myelin
% % %       Xz=ALLCM{jj}{CC}; Xstr=CMLBLnx{jj}{CC}; % Xnx=ALLCMnx{jj}{CC};
% % % 
% % %     % Combine all X predictors
% % %       useD={Xz  Xdz};     DLBL={Xstr  DISTLBL}; % useDnx={Xnx  Xd};
% % %       useN=length(useD); ixd=useN; nd=ixd+1; iy=nd; PLTLBL=DLBL;
% % % 
% % % 
% % %     % Filter for ALL vs CON vs DIS edges
% % %       switch SNP 
% % %           case 0
% % %               Dcdastr=' (ALL-NP)';
% % %           case 1
% % %               for kk=1:ixd;  useD{kk}=useD{kk}(masknz); end
% % %               Dcdastr=' (CON-NP)';
% % %           case 2
% % %               for kk=1:ixd;  useD{kk}=useD{kk}(~masknz); end
% % %               Dcdastr=' (DIS-NP)';
% % %       end
% % % 
% % % 
% % %     % ---- Combine & compute metadata        
% % %     % Combine X & Y
% % %       ALLD=cell(1,nd);   for ii=1:useN; ALLD{ii}=useD{ii};     end; ALLD{iy}=Yz_m; 
% % %       % ALLDnx=cell(1,nd); for ii=1:useN; ALLDnx{ii}=useDnx{ii}; end; ALLDnx{iy}=Y_m;
% % % 
% % %     % ---- MODEL SPECS
% % %     % Model specs 1: FC ~ 1 + CMi_Wj + ED
% % %       mdlstrgen='y~1+x1+x2';
% % %       mdlstrlonggen=mdlstrgen;
% % %       pltlbl_mdl=['FC ~ 1 + ' PLTLBL{1} ' + ' PLTLBL{2}];  
% % %       D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{1};
% % %       Dt=rmzeros(D);  X=Dt(:,1:ixd); y=Dt(:,iy);                              % Remove 0-valued edges
% % %     % Modeling
% % %       mdlspec=[ylbl{1} '~1+' PLTLBL{1} '+' PLTLBL{2}];
% % %       M=fitlm(X,y,mdlspec,'VarNames',[PLTLBL ylbl{1}]);
% % %       PLTLBL_hi{jj}=M.Coefficients.Properties.RowNames(2:end)';
% % % 
% % % 
% % %     % ---- Nodal model
% % %       for ii = 1 : Nfc
% % %           itfail=0; itpass=0; it=0;
% % %         % Extract all X data with the FC data for this round
% % %           D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
% % %         % Cut the subcortical data out of the SC if not present in FC
% % %           if ~all(cellfun(@(c) length(c)==Nnode,D))
% % %               it=Nstx;
% % %               D(1:ixd)=cellfun(@(c) c(it+1:end,it+1:end), D(1:ixd),'UniformOutput',0);
% % %           end
% % %         % Number of nodes for this round
% % %           tNnode = unique(cellfun(@(c) length(c),D));
% % %           if length(tNnode)>1; disp('WARNING: MISMATCHED NODE COUNTS'); keyboard; end
% % %         % Specify model
% % %           mdlspec=[ylbl{ii} '~1+' PLTLBL{1} '+' PLTLBL{2}]; 
% % %           mdlstrfull=mdlspec;
% % %           disp(['Running model for y-data: ' ylbl{ii}])
% % %           for nn = 1 : tNnode
% % %               % Select subset of data (individual nodes)
% % %                 ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(:,nn); end
% % %               % rm zeros
% % %                 Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
% % %               % Ensure N data points > Ncoeffs
% % %                 if size(X,1)<=Ncoeff+10
% % %                     disp(['Node fail: ' num2str(nn)]); itfail=itfail+1;
% % %                     nullcfs1=nan(1,Ncoeff);            nullcfs2=nan(1,Ncoeff);                zerocfs1=zeros(1,Ncoeff);
% % %                     R2_node(nn+it,ii,jj)=nan;          DA_node_prcnt(:,nn+it,ii,jj)=nullcfs2; B_node(:,nn+it,ii,jj)=zerocfs1;
% % %                     PVAL_node(:,nn+it,ii,jj)=nullcfs1; DA_node_value(:,nn+it,ii,jj)=nullcfs2; 
% % %                                                        DA_node_ranks(:,nn+it,ii,jj)=nullcfs2;
% % %                     continue
% % %                 end
% % %                 itpass=itpass+1;
% % %               % Model
% % %                 M=fitlm(X,y,mdlspec,'VarNames',[PLTLBL ylbl{ii}]);
% % %               % Dominance Analysis
% % %                 [tmp,~]=dominance_interaction(X,y,mdlstrgen);
% % %                 DA_node_prcnt(:,nn+it,ii,jj)=tmp.Dominance_percentage; 
% % %                 DA_node_value(:,nn+it,ii,jj)=tmp.Dominance_value;
% % %                 DA_node_ranks(:,nn+it,ii,jj)=tmp.Rank;
% % %               % Save Coeffs
% % %                 B_node(:,nn+it,ii,jj)=M.Coefficients.Estimate(2:end);  R2_node(nn+it,ii,jj)=M.Rsquared.Adjusted; 
% % %                 PVAL_node(:,nn+it,ii,jj)=M.Coefficients.pValue(2:end); DFE_node(nn+it,ii,jj)=M.DFE;
% % %               % Store empirical & predicted y values for each node
% % %                 y_pred_node{nn,ii,jj}=[y M.Fitted];  corr_pred_node(nn,ii,jj)=corr(y,M.Fitted);
% % %           end
% % %           disp([sprintf('Node pass rate = %.1f',(itpass / tNnode * 100)) '%']);
% % %       end
% % % 
% % %     % Pool Pred & Emp y for all nodes
% % %       for ii=1:Nfc
% % %           Dt=y_pred_node(:,ii,jj); Dt(cellfun(@(c) isempty(c),Dt))=[];
% % %           Dx=cell2mat(cellfun(@(c) c(:,1),Dt,'UniformOutput',0));
% % %           Dy=cell2mat(cellfun(@(c) c(:,2),Dt,'UniformOutput',0));
% % %           corr_pred_node_pool(jj,ii)=corr(Dx,Dy,'type','Pearson');
% % %       end
% % %   end
% % %   R2_node(R2_node<0)=0;  R2_node(isnan(R2_node))=0;
% % % 
% % % 
% % % 
% % % %% ------- Across model resolutions
% % % 
% % % % --- Simulated vs Empirical (pooled edges bar plot)
% % %   ppostmp=getPlotPos([-2559 878 1022 330],2,0);
% % %   for jj=1:Nw8s
% % %       D=[corr_pred_g(jj,:)' corr_pred_ntwk_pool(jj,:)' corr_pred_node_pool(jj,:)'];
% % %       myfig([PLTLBL_hi{jj}{1} ' ' str],ppostmp(jj,:)); B=bar(ylbl,D,0.6); 
% % %       ylabel('Pearson''s r'); title('empirical vs predicted'); font(fntsz,usefont); 
% % %       cmap_tmp=gray(3); for ii=1:3; B(ii).FaceColor=cmap_tmp(ii,:); end
% % %       legend({'global' 'network' 'node'}); ylim([0 1]);
% % %       grid on; set(gca,'XGrid','off','YTick',0:.2:1);
% % %   end
% % % 
% % % 
% % % 
% % % 
% % % 
% % % %% =============       Modeling 1b: FCx ~ CMi_Wj           ============== %
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % Modeling is performed at node & network levels only using all CMs 
% % % % (one at a time) computed from a single edge weight.
% % % % The goal is to evaluate the spatial patterns of myelin communication.
% % % 
% % % Ncoeff=1; 
% % % SNP=0;                   % switch for ALL (0) vs CON (1) vs DIS (2) edges
% % % jj=1;     % 1 is caliber, 2 is myelin
% % % 
% % %     % % Euclidean distance
% % %     %   Xd=EDnrm; DISTLBL='ED'; Xdz=nzzscore(Xd,0);
% % % 
% % %  corr_samap_beta=zeros(Ncm,Nfc);          
% % %  corr_samap_R2=zeros(Ncm,Nfc);
% % % 
% % % for CC=1:Ncm
% % % 
% % %     disp('-------------------------------------------------')
% % %     disp([' *-*-*   Running:  ' CMLBLnx{jj}{CC} '   *-*-*'])
% % %     disp('-------------------------------------------------')
% % % 
% % % 
% % %     % NETWORK model variables
% % %       R2_ntwk=zeros(Nntwk,Nntwk,Nfc);          DA_ranks_ntwk=zeros(Nntwk,Nntwk,Ncoeff,Nfc);
% % %       DFE_ntwk=zeros(Nntwk,Nntwk,Nfc);         DA_value_ntwk=zeros(Nntwk,Nntwk,Ncoeff,Nfc);
% % %       RMSE_ntwk=zeros(Nntwk,Nntwk,Nfc);        DA_prcnt_ntwk=zeros(Nntwk,Nntwk,Ncoeff,Nfc);         
% % %       B_ntwk=zeros(Nntwk,Nntwk,Ncoeff,Nfc);    PVAL_ntwk=zeros(Nntwk,Nntwk,Ncoeff,Nfc);  
% % %       y_pred_ntwk=cell(Nntwk,Nntwk,Nfc);       corr_pred_ntwk=zeros(Nntwk,Nntwk,Nfc);
% % %       corr_pred_ntwk_pool=zeros(1,Nfc);
% % % 
% % %     % NODE model variables
% % %       B_node=zeros(Ncoeff,Nnode,Nfc);          R2_node=zeros(Nnode,Nfc); 
% % %       PVAL_node=zeros(Ncoeff,Nnode,Nfc);       DFE_node=zeros(Nnode,Nfc); 
% % %       DA_node_value=zeros(Ncoeff,Nnode,Nfc);   DA_node_prcnt=zeros(Ncoeff,Nnode,Nfc);
% % %       DA_node_ranks=zeros(Ncoeff,Nnode,Nfc);   
% % %       y_pred_node=cell(Nnode,Nfc);             corr_pred_node=zeros(Nnode,Nfc);
% % %       corr_pred_node_pool=zeros(1,Nfc);
% % % 
% % % 
% % % 
% % % 
% % %     % ===== Loop over weights
% % %       % for jj=1:Nw8s
% % % 
% % %         % Communication model for caliber or myelin
% % %           Xz=ALLCM{jj}{CC}; Xnx=ALLCMnx{jj}{CC}; Xstr=CMLBLnx{jj}{CC}; Xstr=strrep(Xstr,'-','');
% % % 
% % %         % Combine all X predictors
% % %           % useD={Xz  Xdz};     useDnx={Xnx  Xd};    DLBL={Xstr  DISTLBL};  % with ED
% % %           useD={Xz};     useDnx={Xnx};    DLBL={Xstr};                      % no ED
% % %           useN=length(useD); ixd=useN; nd=ixd+1; iy=nd; PLTLBL=DLBL;
% % % 
% % % 
% % %         % Filter for ALL vs CON vs DIS edges
% % %           switch SNP 
% % %               case 0
% % %                   Dcdastr=' (ALL-NP)';
% % %               case 1
% % %                   for kk=1:ixd;  useD{kk}=useD{kk}(masknz); end
% % %                   Dcdastr=' (CON-NP)';
% % %               case 2
% % %                   for kk=1:ixd;  useD{kk}=useD{kk}(~masknz); end
% % %                   Dcdastr=' (DIS-NP)';
% % %           end
% % % 
% % % 
% % % 
% % %         % ---- Combine & compute metadata
% % %         % Combine X & Y
% % %           ALLD=cell(1,nd);   for ii=1:useN; ALLD{ii}=useD{ii};     end; ALLD{iy}=Yz_m; 
% % %           ALLDnx=cell(1,nd); for ii=1:useN; ALLDnx{ii}=useDnx{ii}; end; ALLDnx{iy}=Y_m;
% % % 
% % %         % ---- MODEL SPECS
% % %         % Model specs 2: FC ~ 1 + CMi_Wj
% % %           mdlstrgen='y~1+x1';
% % %           mdlstrlonggen=mdlstrgen;
% % %           pltlbl_mdl=['FC ~ 1 + ' PLTLBL{1}];  
% % %           D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{1};
% % %           Dt=rmzeros(D);  X=Dt(:,1:ixd); y=Dt(:,iy);                              % Remove 0-valued edges
% % %         % Modeling
% % %           mdlspec=[ylbl{1} '~1+' PLTLBL{1}];
% % %           M=fitlm(X,y,mdlspec,'VarNames',[PLTLBL ylbl{1}]);
% % %           PLTLBL_hi=M.Coefficients.Properties.RowNames(2:end)';
% % % 
% % %         % % % % Model specs 1: FC ~ 1 + CMi_Wj + ED
% % %         % % %   mdlstrgen='y~1+x1+x2';
% % %         % % %   mdlstrlonggen=mdlstrgen;
% % %         % % %   pltlbl_mdl=['FC ~ 1 + ' PLTLBL{1} ' + ' PLTLBL{2}];  
% % %         % % %   D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{1};
% % %         % % %   Dt=rmzeros(D);  X=Dt(:,1:ixd); y=Dt(:,iy);                              % Remove 0-valued edges
% % %         % % % % Modeling
% % %         % % %   mdlspec=[ylbl{1} '~1+' PLTLBL{1} '+' PLTLBL{2}];
% % %         % % %   M=fitlm(X,y,mdlspec,'VarNames',[PLTLBL ylbl{1}]);
% % %         % % %   PLTLBL_hi{jj}=M.Coefficients.Properties.RowNames(2:end)';
% % % 
% % % 
% % % 
% % %         % ---- Network-level Modeling
% % % 
% % %           for ii = 1 : Nfc
% % %             % Extract all X data with the FC data for this round
% % %               D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
% % %             % Cut the subcortical data out of the SC if not present in FC
% % %               it=0; usepinfo=pinfo;
% % %               if ~all(cellfun(@(c) length(c)==Nnode,D))
% % %                   it=1; usepinfo=pinfo2;
% % %                   D(1:ixd)=cellfun(@(c) c(Nstx+1:end,Nstx+1:end), D(1:ixd),'UniformOutput',0);
% % %               end
% % %               useNntwk=length(usepinfo.clabels);
% % %             % Specify model
% % %               mdlspec=[ylbl{ii} '~1+' PLTLBL{1}]; 
% % %               disp(['Running model for y-data: ' ylbl{ii}])
% % % 
% % %               for iio = 1 : useNntwk
% % %                   nio=usepinfo.cis{iio};
% % %                   for iii = iio : useNntwk
% % %                       nii=usepinfo.cis{iii};
% % %                       ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(nio,nii); end 
% % %                       Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
% % %                       try  
% % %                           M=fitlm(X,y,mdlspec,'VarNames',[PLTLBL ylbl{ii}]); 
% % %                           if M.DFE <= 9; error(['Insufficient DFE: ' usepinfo.clabels{iio} '-' usepinfo.clabels{iii}]); end
% % %                           R2_ntwk(iio+it,iii+it,ii)=M.Rsquared.Adjusted;      PVAL_ntwk(iio+it,iii+it,:,ii)=M.Coefficients.pValue(2:end);
% % %                           RMSE_ntwk(iio+it,iii+it,ii)=M.RMSE;
% % %                           DFE_ntwk(iio+it,iii+it,ii)=M.DFE;                   B_ntwk(iio+it,iii+it,:,ii)=M.Coefficients.Estimate(2:end); 
% % %                           % [tmp,~]=dominance_interaction(X,y,mdlstrgen);              DA_ranks_ntwk(iio+it,iii+it,:,ii)=tmp.Rank;
% % %                           % DA_value_ntwk(iio+it,iii+it,:,ii)=tmp.Dominance_value;     DA_prcnt_ntwk(iio+it,iii+it,:,ii)=tmp.Dominance_percentage; 
% % %                         % Store empirical & predicted y values
% % %                           y_pred_ntwk{iio+it,iii+it,ii}=[y M.Fitted];                   corr_pred_ntwk(iio+it,iii+it,ii)=corr(y,M.Fitted);
% % %                       catch
% % %                           disp(['FAILED BLOCK: Network-pair ' usepinfo.clabels{iio} '-' usepinfo.clabels{iii}])
% % %                           DA_value_ntwk(iio+it,iii+it,:,ii)=zeros(Ncoeff,1);   R2_ntwk(iio+it,iii+it,ii)=0;               
% % %                           DA_ranks_ntwk(iio+it,iii+it,:,ii)=zeros(Ncoeff,1);   B_ntwk(iio+it,iii+it,:,ii)=zeros(Ncoeff,1);              
% % %                           DA_prcnt_ntwk(iio+it,iii+it,:,ii)=zeros(Ncoeff,1);   RMSE_ntwk(iio+it,iii+it,ii)=inf; 
% % %                           DFE_ntwk(iio+it,iii+it,ii)=0;                        PVAL_ntwk(iio+it,iii+it,:,ii)=zeros(Ncoeff,1);
% % % 
% % %                       end
% % %                   end
% % %               end
% % %           end
% % % 
% % %         % % % % Pool Pred & Emp y for all nodes
% % %         % % %   for ii=1:Nfc
% % %         % % %       Dt=y_pred_ntwk(:,:,ii); Dt=Dt(:); Dt(cellfun(@(c) isempty(c),Dt))=[];
% % %         % % %       Dx=cell2mat(cellfun(@(c) c(:,1),Dt,'UniformOutput',0));
% % %         % % %       Dy=cell2mat(cellfun(@(c) c(:,2),Dt,'UniformOutput',0));
% % %         % % %       corr_pred_ntwk_pool(ii)=corr(Dx,Dy,'type','Pearson');
% % %         % % %   end
% % %       % end
% % % 
% % % 
% % %     % ---- Plotting Network-Level (pairwise)
% % %     % Pairwise Networks
% % %       str2=['Network-level' Dcdastr];
% % % 
% % %     % % % % --- rsNetworks on surface
% % %     % % % % Surface plot
% % %     % % %   D=[1:Nntwk]'; useN=1; TMPLBL='RSNs'; TMPSTR=['Yeo 7-rsNetworks: ' str]; usecmap=[ .75 .75 .75 ; cmap_ntwk];
% % %     % % %   D2=nan(Nnode,useN); for ii=1:Nntwk; for jj=1:useN; D2(pinfo.cis{ii},jj)=D(ii,jj); end; end
% % %     % % %   Hcx=plot_conn_surf(D2,pinfo,'cortex',TMPLBL, TMPSTR); setsurf(Hcx,[0 Nntwk],usecmap);
% % %     % % % % Legend
% % %     % % %   myfig; for ii=1:Nntwk; plot([ii ii],1:2,'-','LineWidth',10,'Color',cmap_ntwk(ii,:)); hold on; end
% % %     % % %   legend(LBL_ntwk); font(fntsz,usefont);
% % %     % % % 
% % % 
% % %     % ----------- Matrix Plots
% % %       pposvals=[498  119 -260];
% % %       axopt={'XTick',1:Nntwk,'XTickLabel',LBL_ntwk,'XTickLabelRotation',40, ...
% % %              'YTick',1:Nntwk,'YTickLabel',LBL_ntwk,'YTickLabelRotation',40,'FontSize',fntsz};
% % %       meshopt={'EdgeColor','k','LineWidth',1.5};
% % %     % --- R2
% % %       D=R2_ntwk; D(D<0)=0; TMPSTR='R^2: '; usecmap=cm3; CLIM=[0 ceil(max(D(:))*10)/10];
% % %       % for jj=1:Nw8s
% % %           TMPSTR2=[TMPSTR PLTLBL_hi{1}]; ppostmp2=getPlotPos([-2559 pposvals(1) 315 299],Nfc,0);
% % %           for ii=1:Nfc
% % %             % Matrix with blank upper tri
% % %               D2=D(:,:,ii)'; tmask=tril(true(size(D2))); myfig([TMPSTR2 ', ' ylbl{ii}],ppostmp2(ii,:));
% % %               imagesc(D2,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; clim(CLIM); font(fntsz,usefont);
% % %               set(gca,axopt{:}); font(fntsz,usefont); hold on; % title(ylbl{ii}); 
% % %               for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % %               set(gca,'XTickLabel','','YTickLabel',''); title '';  axis off; %colorbar off;
% % %           end
% % %       % end
% % %     % --- Beta (p<0.05)
% % %       D=B_ntwk; D(PVAL_ntwk>0.05)=0; TMPSTR='beta: '; usecmap=cm2;  
% % %       % for jj=1:Nw8s
% % %           Dt=squeeze(D(:,:,1,:)); CMAX=round(prctile(abs(Dt(:)),99),2); CLIM=[CMAX*-1  CMAX];
% % %           TMPSTR2=[TMPSTR PLTLBL_hi{1}]; ppostmp2=getPlotPos([-2559 pposvals(2) 315 299],Nfc,0);
% % %           for ii=1:Nfc
% % %               D2=Dt(:,:,ii)'; tmask=tril(true(size(D2))); myfig([TMPSTR2 ', '  ylbl{ii}],ppostmp2(ii,:));
% % %               imagesc(D2,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; clim(CLIM); font(fntsz,usefont);
% % %               set(gca,axopt{:}); font(fntsz,usefont); hold on; % title(ylbl{ii}); 
% % %               for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % %               set(gca,'XTickLabel','','YTickLabel',''); title ''; axis off; % colorbar off;
% % %           end
% % %       % end
% % %     % Beta (p<0.05) MEDIAN across FC 
% % %       D=B_ntwk; D(PVAL_ntwk>0.05)=0; TMPSTR='Median \beta_p_<_0_._0_5: '; usecmap=cm2;
% % %       ppostmp=getPlotPos([-2559 -312 375 351],2);
% % %       % for jj=1:Nw8s 
% % %           Dt=squeeze(D(:,:,1,:)); TMPSTR2=[TMPSTR PLTLBL_hi{1}]; 
% % %           Dt(Dt==0)=nan; Dt=nanmedian(Dt,3); Dt(isnan(Dt))=0; Dt=Dt'; %CLIM=[-1 1];
% % %           CMAX=round(prctile(abs(Dt(:)),99),2); CLIM=[CMAX*-1  CMAX];
% % %           tmask=tril(true(size(Dt))); myfig(TMPSTR2,ppostmp(1,:)); 
% % %           imagesc(Dt,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; font(fntsz,usefont); 
% % %           set(gca,axopt{:}); font(fntsz,usefont); hold on; clim(CLIM);
% % %           for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % %           set(gca,'XTickLabel','','YTickLabel',''); title ''; axis off; %colorbar off;
% % %       % end
% % % 
% % %     % --- Dominance rankings
% % %     % % %   D=DA_ranks_ntwk; D(D==0)=nan; D=(D.*-1)+Ncoeff+1; D(isnan(D))=0; 
% % %     % % %   TMPSTR='Dominance Rank: '; usecmap=cm4;  CLIM=[0 Ncoeff];
% % %     % % %   % for jj=1:Nw8s
% % %     % % %       Dt=squeeze(D(:,:,1,:)); TMPSTR2=[TMPSTR PLTLBL_hi{1}]; ppostmp2=getPlotPos([-2559 pposvals(3) 315 299],Nfc,0); 
% % %     % % %       for ii=1:Nfc
% % %     % % %           D2=Dt(:,:,ii)'; tmask=tril(true(size(D2))); myfig([TMPSTR2 ', ' ylbl{ii}],ppostmp2(ii,:));
% % %     % % %           imagesc(D2,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; clim(CLIM); font(fntsz,usefont);
% % %     % % %           set(gca,axopt{:}); font(fntsz,usefont); hold on; % title(ylbl{ii}); 
% % %     % % %           for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % %     % % %           set(gca,'XTickLabel','','YTickLabel',''); title ''; axis off; %colorbar off;
% % %     % % %       end
% % %     % % %   % end
% % %     % % % % Dominance rank MEDIAN across FC
% % %     % % %   D=DA_ranks_ntwk; D(D==0)=nan; D=(D.*-1)+Ncoeff+1; usecmap=cm4; CLIM=[0 Ncoeff];
% % %     % % %   ppostmp=getPlotPos([-2559 -312 375 351],2); TMPSTR='Median dominance rank: ';
% % %     % % %   % for jj=1:Nw8s 
% % %     % % %       D2=squeeze(D(:,:,1,:));  TMPSTR2=[TMPSTR PLTLBL_hi{1}];
% % %     % % %       DPLT=nanmedian(D2,3); DPLT(isnan(DPLT))=0; DPLT=DPLT'; 
% % %     % % %       tmask=tril(true(size(DPLT))); myfig(TMPSTR2,ppostmp(1,:));
% % %     % % %       imagesc(DPLT,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; clim(CLIM); font(fntsz,usefont);
% % %     % % %       set(gca,axopt{:}); font(fntsz,usefont); hold on;
% % %     % % %       for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % %     % % %       set(gca,'XTickLabel','','YTickLabel',''); title ''; colorbar off; axis off;
% % %     % % %   end
% % %     % % % 
% % % 
% % %     % % % % ----- Empirical vs Predicted:
% % %     % % % % Scatter plot: all edges pooled
% % %     % % %   % % % pposvals=[498  119];
% % %     % % %   % % % for jj=1:Nw8s
% % %     % % %   % % %     ppostmp2=getPlotPos([-2559 pposvals(jj) 315 299],Nfc,0);
% % %     % % %   % % %     for ii=1:Nfc
% % %     % % %   % % %         Dt=y_pred_ntwk(:,:,ii,jj); Dt=Dt(:); Dt(cellfun(@(c) isempty(c),Dt))=[]; 
% % %     % % %   % % %         Dx=cell2mat(cellfun(@(c) c(:,1),Dt,'UniformOutput',0));
% % %     % % %   % % %         Dy=cell2mat(cellfun(@(c) c(:,2),Dt,'UniformOutput',0));
% % %     % % %   % % %         myfig([PLTLBL_hi{jj}{1} str2],ppostmp2(ii,:)); heatscatter(Dx,Dy); colormap(gray);
% % %     % % %   % % %         xlabel('empirical'); ylabel('predicted'); font(fntsz,usefont); 
% % %     % % %   % % %         title([sprintf('r=%.2f (',corr_pred_ntwk_pool(jj,ii)) ylbl{ii} ')']);
% % %     % % %   % % %         colorbar off;
% % %     % % %   % % %     end
% % %     % % %   % % % end
% % %     % % % % Bar plot: edges pooled
% % %     % % %   ppostmp=getPlotPos([-2559 467 467 330],2,0);
% % %     % % %   for jj=1:Nw8s
% % %     % % %       myfig([PLTLBL_hi{jj}{1} ', ' str2],ppostmp(jj,:)); B=bar(ylbl,corr_pred_ntwk_pool(jj,:),0.5); 
% % %     % % %       ylabel('Pearson''s r'); title('empirical vs predicted'); font(fntsz,usefont); grid on;
% % %     % % %       B.CData=cmap_pred(jj,:); B.FaceColor='flat'; ylim([0 1]);
% % %     % % %   end
% % % 
% % %  pause
% % %  close all
% % % 
% % % 
% % %     % -------------------------------------------------------------------------
% % %     %% =======================        Node-level       ========================
% % %     % -------------------------------------------------------------------------
% % % 
% % % 
% % %         % ---- Nodal model
% % %           for ii = 1 : Nfc
% % %               itfail=0; itpass=0; it=0;
% % %             % Extract all X data with the FC data for this round
% % %               D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
% % %             % Cut the subcortical data out of the SC if not present in FC
% % %               if ~all(cellfun(@(c) length(c)==Nnode,D))
% % %                   it=Nstx;
% % %                   D(1:ixd)=cellfun(@(c) c(it+1:end,it+1:end), D(1:ixd),'UniformOutput',0);
% % %               end
% % %             % Number of nodes for this round
% % %               tNnode = unique(cellfun(@(c) length(c),D));
% % %               if length(tNnode)>1; disp('WARNING: MISMATCHED NODE COUNTS'); keyboard; end
% % %             % Specify model
% % %               mdlspec=[ylbl{ii} '~1+' PLTLBL{1}]; 
% % %               mdlstrfull=mdlspec;
% % %               disp(['Running model for y-data: ' ylbl{ii}])
% % %               for nn = 1 : tNnode
% % %                   % Select subset of data (individual nodes)
% % %                     ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(:,nn); end
% % %                   % rm zeros
% % %                     Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
% % %                   % Ensure N data points > Ncoeffs
% % %                     if size(X,1)<=Ncoeff+10
% % %                         disp(['Node fail: ' num2str(nn)]); itfail=itfail+1;
% % %                         nullcfs1=nan(1,Ncoeff);            nullcfs2=nan(1,Ncoeff);                zerocfs1=zeros(1,Ncoeff);
% % %                         R2_node(nn+it,ii)=nan;             DA_node_prcnt(:,nn+it,ii)=nullcfs2;    B_node(:,nn+it,ii)=zerocfs1;
% % %                         PVAL_node(:,nn+it,ii)=nullcfs1;    DA_node_value(:,nn+it,ii)=nullcfs2; 
% % %                                                            DA_node_ranks(:,nn+it,ii)=nullcfs2;
% % %                         continue
% % %                     end
% % %                     itpass=itpass+1;
% % %                   % Model
% % %                     M=fitlm(X,y,mdlspec,'VarNames',[PLTLBL ylbl{ii}]);
% % %                   % % % % Dominance Analysis
% % %                   % % %   [tmp,~]=dominance_interaction(X,y,mdlstrgen);
% % %                   % % %   DA_node_prcnt(:,nn+it,ii)=tmp.Dominance_percentage; 
% % %                   % % %   DA_node_value(:,nn+it,ii)=tmp.Dominance_value;
% % %                   % % %   DA_node_ranks(:,nn+it,ii)=tmp.Rank;
% % %                   % Save Coeffs
% % %                     B_node(:,nn+it,ii)=M.Coefficients.Estimate(2:end);  R2_node(nn+it,ii)=M.Rsquared.Adjusted; 
% % %                     PVAL_node(:,nn+it,ii)=M.Coefficients.pValue(2:end); DFE_node(nn+it,ii)=M.DFE;
% % %                   % Store empirical & predicted y values for each node
% % %                     % y_pred_node{nn,ii}=[y M.Fitted];  corr_pred_node(nn,ii)=corr(y,M.Fitted);
% % %               end
% % %               disp([sprintf('Node pass rate = %.1f',(itpass / tNnode * 100)) '%']);
% % %           end
% % % 
% % %         % Pool Pred & Emp y for all nodes
% % %           % % % for ii=1:Nfc
% % %           % % %     Dt=y_pred_node(:,ii); Dt(cellfun(@(c) isempty(c),Dt))=[];
% % %           % % %     Dx=cell2mat(cellfun(@(c) c(:,1),Dt,'UniformOutput',0));
% % %           % % %     Dy=cell2mat(cellfun(@(c) c(:,2),Dt,'UniformOutput',0));
% % %           % % %     corr_pred_node_pool(ii)=corr(Dx,Dy,'type','Pearson');
% % %           % % % end
% % %       % end
% % %       R2_node(R2_node<0)=0;  R2_node(isnan(R2_node))=0;
% % % 
% % %     % =======  Plot node results
% % % 
% % %     % ----- R2: 
% % %     % Data for CTX
% % %       D2=R2_node; D2(D2<.01)=.01; TMPSTR=['R^2: ' PLTLBL_hi{1} str]; CLIM=[0 ceil(max(D2(:))*10)/10]; usecmap=[.75 .75 .75; cm4]; 
% % %     % CTX surfaces
% % %       Hcx=plot_conn_surf(D2,pinfo2,'cortex',ylbl,TMPSTR); setsurf(Hcx,CLIM,usecmap);
% % %       Hcx{1}.figure.Position=[-1.7778 -0.0156 0.9000 0.9000];
% % %       Hcx{2}.figure.Position=[-1.1306 -0.0156 0.9000 0.9000];
% % % 
% % %     % --- Correlation of R2 with S-A map
% % %       D2=R2_node; dx=SA_map; dxlbl='S-A map'; savecorr=zeros(1,Nfc); % opt={pltlbl_mdl,cm};% heatscat settings
% % %       ppostmp=getPlotPos([-2559 -396 224 303],Nfc,1); CLIM=[0 ceil(max(D2(:))*10)/10];
% % %     % Scatter plot
% % %       for xx=1:Nfc
% % %           dy=D2(:,xx); dylbl=['R^2 (' ylbl{xx} ')'];
% % %           ddisp_corr(dx,dy,dxlbl,dylbl,ylbl{xx},cm,ppostmp(xx,:),fntsz); 
% % %           set(gca,'XTickLabels',''); colorbar off; ylim(CLIM); ylabel '';
% % %           line(xlim,[0 0],'Color','r','LineStyle','--','LineWidth',2);
% % %         % Save correlation
% % %           savecorr(xx)=corr(dx(dy~=0),dy(dy~=0));
% % %       end
% % %     % Bar plot
% % %       myfig(str2,[-759 -396 404 330]); B=bar(ylbl,savecorr,0.5); % [-1115 359 467 330]
% % %       ylabel('Pearson''s r'); title([ 'R^2 vs ' dxlbl]); font(fntsz,usefont); 
% % %       grid on; set(gca,'XGrid','off','XMinorGrid','off');
% % %       B.CData=CLRS.carolinablue; B.FaceColor='flat'; % ylim([-.52 .02]);
% % %     % Store 
% % %       corr_samap_R2(CC,:)=savecorr;
% % % 
% % % 
% % % pause
% % % close all
% % % 
% % % 
% % %     % ----- BETA values: SC, Myelin & length
% % %     % -- CTX
% % %     % Data
% % %       D2=B_node; D2(PVAL_node > 0.05)=0; usecmap=cm2;
% % %       Dt=squeeze(D2(1,:,:)); CMAX=round(prctile(abs(Dt(Dt~=0)),99),1); CLIM1=[CMAX*-1 CMAX]; TMPSTR1=['Betas: ' PLTLBL_hi{1} str];
% % %       Hcx=plot_conn_surf(Dt,pinfo2,'cortex',ylbl,TMPSTR1); setsurf(Hcx,CLIM1,usecmap);
% % %       Hcx{1}.figure.Position=[-1.7778 -0.0156 0.9000 0.9000];
% % %       Hcx{2}.figure.Position=[-1.1306 -0.0156 0.9000 0.9000];
% % % 
% % %     % --- Correlation of BETAS with S-A map
% % %       D2=B_node; D2(PVAL_node > 0.05)=0; Dt=squeeze(D2(1,:,:)); CLIM1=[floor(min(Dt(Dt~=0))*10)/10  ceil(max(Dt(Dt~=0))*10)/10];
% % %       dx=SA_map; dxlbl='S-A map'; savecorr=zeros(1,Nfc); % opt={pltlbl_mdl,cm};% heatscat settings
% % %       ppostmp=getPlotPos([-2559 -396 224 303],Nfc,1);
% % %     % Scatter plot
% % %       for xx=1:Nfc
% % %           dy=Dt(:,xx); dylbl=['\beta (' ylbl{xx} ')'];
% % %           ddisp_corr(dx,dy,dxlbl,dylbl,['betas, ' ylbl{xx}],cm,ppostmp(xx,:),fntsz); 
% % %           set(gca,'XTickLabels',''); colorbar off; ylim(CLIM1); ylabel '';
% % %           line(xlim,[0 0],'Color','r','LineStyle','--','LineWidth',2);
% % %         % Save correlation
% % %           savecorr(xx)=corr(dx(dy~=0),dy(dy~=0));
% % %       end
% % %     % Bar plot
% % %       myfig(str2,[-759 -396 404 330]); B=bar(ylbl,savecorr,0.5); % [-1115 359 467 330]
% % %       ylabel('Pearson''s r'); title([ '\betas vs ' dxlbl]); font(fntsz,usefont); 
% % %       grid on; set(gca,'XGrid','off','XMinorGrid','off');
% % %       B.CData=CLRS.carolinablue; B.FaceColor='flat'; % ylim([-.52 .02]);
% % %     % Store 
% % %       corr_samap_beta(CC,:)=savecorr;
% % % 
% % %  pause
% % %  close all
% % % 
% % % 
% % % 
% % % end
% % % 
% % % % ---- Plot all correlations with SA map
% % % % --- Matrices
% % %   ppostmp=getPlotPos([-2559 389 474 408],2,0);
% % % % - R2 vs SA
% % %   D=corr_samap_R2; TMPSTR1='R^2: '; TMPSTR2=[TMPSTR1 useDstr{jj}];
% % %   CLIM=[-1 1];% CMAX=ceil(max(abs(D(:)))*10)/10; CLIM=[CMAX*-1 CMAX];
% % %   myfig([TMPSTR2 str],ppostmp(1,:)); imagesc(D); colorbar; colormap(cm2); 
% % %   clim(CLIM); font(fntsz,usefont); % title(TMPSTR2); 
% % %   set(gca,'YTickLabels',CMLBL_simple,'XTick',1:Nfc,'XTickLabels',ylbl,'XTickLabelRotation',45)
% % % % - BETA vs SA
% % %   D=corr_samap_beta; TMPSTR1='betas: '; TMPSTR2=[TMPSTR1 useDstr{jj}];
% % %   CLIM=[-1 1];% CMAX=ceil(max(abs(D(:)))*10)/10; CLIM=[CMAX*-1 CMAX];
% % %   myfig([TMPSTR2 str],ppostmp(2,:)); imagesc(D); colorbar; colormap(cm2); 
% % %   clim(CLIM); font(fntsz,usefont); % title(TMPSTR2); 
% % %   set(gca,'YTickLabels',CMLBL_simple,'XTick',1:Nfc,'XTickLabels',ylbl,'XTickLabelRotation',45)





