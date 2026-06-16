% Exploratory analysis for paper 3
%
% Final hierarchical regression analysis: adding caliber to base model
%
% 1. Modeling 1a:   Hiearchical
%       Base model: FCx ~ CMi_binary + CMi_caliber + ED
%       Test predictors: MTsat, g-ratio & delay
%       Level 1: individual predictors
%       Level 2: pairwise combinations of predictors
%       Level 3: all predictors
%       + Results confirmed using residualization
%       + all edges together
%
% 2. Modeling 1b: Blockwise RSN prediction (UNFINISHED)
% 3. Modeling 1c: Dual BICS
% 4. Modeling 1d: Dual BICS ΔΔR² (myelin vs caliber)
%
% 2025 Mark C Nelson MNI
%--------------------------------------------------------------------------

%% Params
% Data
  MCstr         = 'MTsat-tm';                                               % {'MTsat' 'R1-tm' 'gratio-tm' 'gratio-ts'}
  MC2str        = 'gratio-ts';                                               % {'MTsat' 'R1-tm' 'gratio-tm' 'gratio-ts'}
  % SCstr         = 'COMMITscl';                                              % {'COMMITscl' 'COMMITsclvol' 'COMMITsclNoS'};
  SCstr         = 'COMMITsclNoS';                                              % {'COMMITscl' 'COMMITsclvol' 'COMMITsclNoS'};
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
  % cmap_pred=[CLRS.carrot; CLRS.aztcprpl; .8,.8,.8;  .50,.50,.50];              % Predictors; COMMIT=red; MYLNPC1=blue; LoS=lightGray; ED=darkGray
  % cmap_pred=[CLRS.carrot; CLRS.aztcprpl; CLRS.livrplgreen; .8,.8,.8];            % Predictors
  cmap_pred=[CLRS.aztcprpl; CLRS.livrplred; CLRS.livrplgreen; .8,.8,.8];            % Predictors
  cmap_5cm =[.24,.74,.71; .64,.40,.84; .55,.80,.99; .95,.82,.18; .98,.51,.45]; % 5-CM; SPE=green; NE=purple; SI=blue; CMY=gold DE=red; (close to Seguin2020)
  cmap_6cm =[.24,.74,.71; .64,.40,.84; .55,.80,.99; .59,.47,.45; .95,.82,.18; .98,.51,.45]; % 6-CM; SPE=green; NE=purple; SI=blue; PT=brown; CMY=gold DE=red; (close to Seguin2020)
  cmap_intract=[.13,.67,.58;   .80,.12,.25];
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

% Yeo 7-Networks 
  Dt=[1:Nntwk]'; D2=nan(Nnode,1); for ii=1:Nntwk; D2(pinfo.cis{ii},1)=Dt(ii,1); end; YeoRSNs=D2;


%% ================  Load & prep EDGE data 

  loaddir_edge           = [locDeriv '/0_finalData_edges/grp_' parc '_CBthr-' num2str(thr) '_distDep/'];
  
% --- Load EDGE data
  load([loaddir_edge SCstr '.mat'],'Dtg');                  SC=Dtg;         % Connection-Strength
  load([loaddir_edge MCstr '.mat'],'Dtg');                  MC=Dtg;         % Myelin
  load([loaddir_edge MC2str '.mat'],'Dtg');                 MC2=Dtg;        % Myelin 2
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


% --- masks for CON vs DIS edges
  masknz    = SC~=0;                                                        % mask for all non zero edges in connectomes (binary connectome NOT uniform across weights)
  maskdis=~masknz; maskdis(1:Nnode+1:end)=0;                                % mask for edges with SC=0

% ---- Normalize
  SCnrm     = SC   ./ ( max(SC(:))  + (max(SC(:))*.01)  );                  % to avoid 0 lengths at existing edges
  MCnrm     = MC   ./ ( max(MC(:))  + (max(MC(:))*.01)  );
  MC2nrm    = MC2  ./ ( max(MC2(:)) + (max(MC2(:))*.01) );
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
%% ============      Communication Models (EDGE data)       ============= %
% -------------------------------------------------------------------------


% Setup
  str1='Communication Models: ';

% ---- Define the data
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




% ==========   Visualize CMs general trends
% Caliber vs myelin for each CM
% vs BOLD-FC
% vs ED

% % %   xx=3; 
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



% ==========   CMs vs all FC

  % % % xx=6; 
% % %   ppostmp=[-539 -349 -160 28 218 407 597]; useY=Yz_m([8:-1:3 1]); YLBL=ylbl([8:-1:3 1]);
% % %   opt4={'XTickLabel',{'all' 'dis' 'con'},'XTickLabelRotation',0}; opt={'XTickLabels','','YTickLabels',''};
% % % 
% % % % ---- zscored data
% % %   D1=ALLCM{1}{xx}; Dstr1=CMLBL{1}{xx}; CMAX=prctile(abs(D1(D1~=0)),99); CLIM1=[CMAX*-1 CMAX];
% % %   D2=ALLCM{2}{xx}; Dstr2=CMLBL{2}{xx}]; CMAX=prctile(abs(D2(D2~=0)),99); CLIM2=[CMAX*-1 CMAX];
% % % % --- Loop over FC
% % %   for yy=1:7
% % %       D3=useY{yy}; Dstr3=YLBL{yy}; CMAX=prctile(abs(D3(D3~=0)),99); CLIM3=[CMAX*-1 CMAX];
% % %       useppos=getPlotPos([-2559 ppostmp(yy) 238 187],10);
% % %     % -- Correlations
% % %       svcorr=zeros(1,6); TMPLBL2=[CMLBL_simple{xx} ' vs ' Dstr3];
% % %       X=D1;          Y=D3;          ii=1; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (ALL)'],[Dstr3 ' (ALL)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end
% % %       X=D2;          Y=D3;          ii=2; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (ALL)'],[Dstr3 ' (ALL)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end; ylabel '';
% % %       X=D1(~masknz); Y=D3(~masknz); ii=3; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (DIS)'],[Dstr3 ' (DIS)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end; ylabel '';
% % %       X=D2(~masknz); Y=D3(~masknz); ii=4; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (DIS)'],[Dstr3 ' (DIS)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end; ylabel '';
% % %       X=D1(masknz);  Y=D3(masknz);  ii=5; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr1 ' (CON)'],[Dstr3 ' (CON)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end; ylabel '';
% % %       X=D2(masknz);  Y=D3(masknz);  ii=6; [~,svcorr(ii)]=ddisp_corr(X,Y,[Dstr2 ' (CON)'],[Dstr3 ' (CON)'],str,cm,useppos(ii,:),20); colorbar off; if yy>1; xlabel ''; end; ylabel '';
% % %     % Summary bar plot
% % %       myfig([Dstr3 str],useppos(7,:)); B=bar(svcorr,0.6);  hold on; % title(TMPLBL2);
% % %       ylabel('\rho'); font(fntsz,usefont); grid on; set(gca,'XGrid','off');
% % %       for pp=1:2:6; B.CData(pp,:)=cmap_pred(1,:); end 
% % %       for pp=2:2:6; B.CData(pp,:)=cmap_pred(2,:); end
% % %       B.FaceColor='flat'; set(gca,'XTick',1.5:2:5.5,opt4{:}); 
% % %       for pp=2:2:4; line([pp+.5 pp+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
% % %     % Connectomes: Yeo RSN ordering
% % %       connplot(D1,[TMPLBL2 str],pinfo,useppos(8,:));  title(Dstr1); clim(CLIM1); colormap(cm2); set(gca,opt{:}); font(18,usefont)
% % %       connplot(D2,[TMPLBL2 str],pinfo,useppos(9,:));  title(Dstr2); clim(CLIM2); colormap(cm2); set(gca,opt{:}); font(18,usefont)
% % %       connplot(D3,[TMPLBL2 str],pinfo,useppos(10,:)); title(Dstr3); clim(CLIM3); colormap(cm2); set(gca,opt{:}); font(18,usefont)
% % %   end


% ==== Box plots of within & between RSNs for CMi caliber vs myelin
% myelin supports greater SPE & NE: between all RSNs, within transmodal
% is myelin-SI also lower in transmodal & between RSNs?
% does PT fit?
% CMY?
% myelin-DE higher in transmodal & between as well?





% ================ SA-axis vs CM-strength
% scatters & bars

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


% -------------------------------------------------------------------------



% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% =======   Modeling 1a: Hiearchical: all edges together
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% + Hiearchical regression
% + base model with BINARY & CALIBER
% + confirmation by residualization


for CC=[6]                              % Switch for communication model
close all; 

    % ----- Setup
      Nit=1;                            % dummy var as dis & con removed here
      Nlvl1=3;                          % number of test predictors
      Nlvl2=Nlvl1*(Nlvl1-1)/2;          % pairwise combinations
      Nlvl3=1;                          % last level
      Ntot=Nlvl1+Nlvl2+Nlvl3;
      useCMLBL=[CMLBL_simple{CC}];
    
    % Euclidean distance
      Xd=EDnrm; DISTLBL='ED'; Xdz=nzzscore(Xd,0);
        
    % --- Communication models
    % myelin & DELAY
      X1z=ALLCM{2}{CC};    X1_str=CMLBLnx{2}{CC};    X1_str=extractBefore(X1_str,'-'); %X1_str=strrep(X1_str,'-','');    % MTsat
      X2z=ALLCM{3}{CC};    X2_str=CMLBLnx{3}{CC};    X2_str=extractBefore(X2_str,'-'); %X2_str=strrep(X2_str,'-','');    % g-ratio
      X3z=ALLCM_delay{CC}; X3_str=CMLBLnx_delay{CC}; X3_str=extractBefore(X3_str,'-'); %X3_str=strrep(X3_str,'-','');    % delay
    
    % Communication model for binary
      Xbz=ALLCM_bin{CC};   Xb_str=CMLBLnx_bin{CC};   Xb_str=strrep(Xb_str,'-','');      % binary
      Xcz=ALLCM{1}{CC};    Xc_str=CMLBLnx{1}{CC};    Xc_str=extractBefore(Xc_str,'-');  % caliber
    
    % % % % TMP PLOTS
    % % %   CLIM=[-3 3];
    % % %   connplot(X1z,X1_str,pinfo); clim(CLIM); colormap(cm2); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk)
    % % %   connplot(X2z,X2_str,pinfo); clim(CLIM); colormap(cm2); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk)
    % % %   connplot(X3z,X3_str,pinfo); clim(CLIM); colormap(cm2); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk)
    % % %   connplot(Xbz,Xb_str,pinfo); clim(CLIM); colormap(cm2); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk)
    % % %   connplot(Xdz,DISTLBL,pinfo); clim(CLIM); colormap(cm2); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk)
    
    
    % Global model variables
      deltaR2_g_ALL=zeros(Nit*Ntot,Nfc);
      pval_g_ALL=zeros(Nit*Ntot,Nfc);
      % % % R2_g_rsd_1_ALL=zeros(Nit*Nlvl1,Nfc);
      % % % R2_g_rsd_2_ALL=zeros(Nit*Nlvl1,Nfc);
      % % % R2_g_rsd_3_ALL=zeros(Nit*Nlvl1,Nfc);
    
    % Network model variables
      deltaR2_ntwk_ALL=zeros(Nit*Ntot,Nntwk,Nntwk,Nfc);
      pval_ntwk_ALL=zeros(Nit*Ntot,Nntwk,Nntwk,Nfc);
      % % % R2_ntwk_rsd_1_ALL=zeros(Nit*Nlvl1,Nntwk,Nntwk,Nfc);
      % % % R2_ntwk_rsd_2_ALL=zeros(Nit*Nlvl1,Nntwk,Nntwk,Nfc);
      % % % R2_ntwk_rsd_3_ALL=zeros(Nit*Nlvl1,Nntwk,Nntwk,Nfc);
    
    % Node model variables
      deltaR2_node_ALL=zeros(Nit*Ntot,Nnode,Nfc);
      pval_node_ALL=zeros(Nit*Ntot,Nnode,Nfc);
    
    % Other
      PLTLBL_ALL=cell(Nit*Ntot,1);
      Dcdastr_ALL=cell(1,Nit);
    
      it=0;
    
    %% LOOP OVER CON & DIS (not actually used)
    for ss=1:Nit
    
        % Combine all X predictors
          useD_base={Xbz     Xcz     Xdz};       useD_test={X1z     X2z     X3z};
          DLBL_base={Xb_str  Xc_str  DISTLBL};   DLBL_test={X1_str  X2_str  X3_str};
          useN_b=length(useD_base);         useN_t=length(useD_test); 
    
        % Filter for ALL vs CON vs DIS edges
          switch ss 
              case 1
                  Dcdastr=' (ALL)';
              % % % case 2 
              % % %     Dcdastr=' (CON)'; 
              % % %     for kk=1:useN_b; useD_base{kk}(~masknz)=0; end
              % % %     for kk=1:useN_t; useD_test{kk}(~masknz)=0; end
              % % % case 3 
              % % %     Dcdastr=' (DIS)'; 
              % % %     for kk=1:useN_b; useD_base{kk}(~maskdis)=0; end
              % % %     for kk=1:useN_t; useD_test{kk}(~maskdis)=0; end
          end    
    
    %% ======== LEVEL 1: each individual edge weight
      
      ixd=useN_b + 1;  iy=ixd+1; 
    
      % ----- LOOP OVER TEST PREDICTORS
    
        for tt = 1 : Nlvl1
    
            it=it+1;
    
            % ---- Combine & compute metadata
            % Combine X & Y
              ALLD=cell(1,iy);            ALLD(1:useN_b)=useD_base;
              ALLD{ixd}=useD_test{tt};    ALLD{iy}=Yz_m; 
              PLTLBL_ALL{it}=[DLBL_test{tt} Dcdastr];
              disp([' --> Running test for: ' PLTLBL_ALL{it}])
    
                
            % ========== Global
    
              pval_g=zeros(1,Nfc);
              deltaR2_g=zeros(1,Nfc);
              % % % R2_g_rsd_1=zeros(1,Nfc);
              % % % R2_g_rsd_2=zeros(1,Nfc);
              % % % R2_g_rsd_3=zeros(1,Nfc);
              for ii = 1 : Nfc
                  disp([' Global: ' ylbl{ii}])
                % Extract all X data with the FC data for this round
                  D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
                % Remove 0-valued edges
                  Dt=rmzeros(D); X=Dt(:,1:ixd); y=Dt(:,iy);
                % Modeling: Hiearchical
                  X_base=X(:,1:useN_b);
                  [deltaR2_g(ii),pval_g(ii)]=compare_models(y,X_base,X);
                % % % % --- Orthogonalization tests
                % % %   X_test=X(:,useN_b+1);
                % % % % 1. Regress binary & ED out of test predictor & rerun model
                % % %   M=fitlm(X_base,X_test);
                % % %   X_resid=M.Residuals.Raw;
                % % %   M=fitlm(X_resid,y);
                % % %   R2_g_rsd_1(ii)=M.Rsquared.Adjusted;
                % % % % 2. Regress binary out of test predictor & rerun model with ED
                % % %   M=fitlm(X_base(:,1),X_test);
                % % %   X_resid=M.Residuals.Raw;
                % % %   M=fitlm([X_resid X_base(:,2)],y);
                % % %   R2_g_rsd_2(ii)=M.Rsquared.Adjusted;
                % % % % 3. Regress binary out of test predictor & rerun model without ED
                % % %   M=fitlm(X_base(:,1),X_test);
                % % %   X_resid=M.Residuals.Raw;
                % % %   M=fitlm(X_resid,y);
                % % %   R2_g_rsd_3(ii)=M.Rsquared.Adjusted;
              end
            % Store outputs
              deltaR2_g_ALL(it,:)=deltaR2_g; 
              pval_g_ALL(it,:)=pval_g;
              % % % R2_g_rsd_1_ALL(it,:)=R2_g_rsd_1;
              % % % R2_g_rsd_2_ALL(it,:)=R2_g_rsd_2;
              % % % R2_g_rsd_3_ALL(it,:)=R2_g_rsd_3;
            
                
        
            % ========== Pairwise networks  
    
              pval_ntwk=zeros(Nntwk,Nntwk,Nfc);
              deltaR2_ntwk=zeros(Nntwk,Nntwk,Nfc);
              % % % R2_ntwk_rsd_1=zeros(Nntwk,Nntwk,Nfc);
              % % % R2_ntwk_rsd_2=zeros(Nntwk,Nntwk,Nfc);
              % % % R2_ntwk_rsd_3=zeros(Nntwk,Nntwk,Nfc);
              for ii = 1 : Nfc
                  disp([' Network: ' ylbl{ii}])
                % Extract all X data with the FC data for this round
                  D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
                  for iio = 1 : Nntwk
                      nio=pinfo.cis{iio};
                      for iii = iio : Nntwk
                          nii=pinfo.cis{iii};
                          ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(nio,nii); end 
                          Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
                        % Modeling: Hiearchical
                          X_base=X(:,1:useN_b);
                          [deltaR2_ntwk(iio,iii,ii),pval_ntwk(iio,iii,ii)]=compare_models(y,X_base,X);
                        % % % % --- Orthogonalization tests
                        % % %   X_test=X(:,useN_b+1);
                        % % % % 1. Regress binary & ED out of test predictor & rerun model
                        % % %   M=fitlm(X_base,X_test);
                        % % %   X_resid=M.Residuals.Raw;
                        % % %   M=fitlm(X_resid,y);
                        % % %   R2_ntwk_rsd_1(iio,iii,ii)=M.Rsquared.Adjusted;
                        % % % % 2. Regress binary out of test predictor & rerun model with ED
                        % % %   M=fitlm(X_base(:,1),X_test);
                        % % %   X_resid=M.Residuals.Raw;
                        % % %   M=fitlm([X_resid X_base(:,2)],y);
                        % % %   R2_ntwk_rsd_2(iio,iii,ii)=M.Rsquared.Adjusted;
                        % % % % 3. Regress binary out of test predictor & rerun model without ED
                        % % %   M=fitlm(X_base(:,1),X_test);
                        % % %   X_resid=M.Residuals.Raw;
                        % % %   M=fitlm(X_resid,y);
                        % % %   R2_ntwk_rsd_3(iio,iii,ii)=M.Rsquared.Adjusted;
    
                      end
                  end
              end
            % Store outputs
              deltaR2_ntwk_ALL(it,:,:,:)=deltaR2_ntwk; 
              pval_ntwk_ALL(it,:,:,:)=pval_ntwk;
              % % % R2_ntwk_rsd_1_ALL(it,:,:,:)=R2_ntwk_rsd_1;
              % % % R2_ntwk_rsd_2_ALL(it,:,:,:)=R2_ntwk_rsd_2;
              % % % R2_ntwk_rsd_3_ALL(it,:,:,:)=R2_ntwk_rsd_3;
    
            
            % ==========        Node-level    
          
              pval_node=zeros(Nnode,Nfc);
              deltaR2_node=zeros(Nnode,Nfc);
              for ii = 1 : Nfc
                  disp([' Node: ' ylbl{ii}])
                % Extract all X data with the FC data for this round
                  D=cell(1,iy); D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
                  for nn = 1 : Nnode
                      % Select subset of data (individual nodes)
                        ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(:,nn); end
                      % rm zeros
                        Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
                      % Modeling: Hiearchical
                        X_base=X(:,1:useN_b);
                        [deltaR2_node(nn,ii),pval_node(nn,ii)]=compare_models(y,X_base,X);
                  end
              end
            % Store outputs
              deltaR2_node_ALL(it,:,:)=deltaR2_node; 
              pval_node_ALL(it,:,:)=pval_node;
    
        end
    
      %% ======== LEVEL 2: pairwise combinations of edge weights
      
        ixd=useN_b + 2;  iy=ixd+1;
        useD_test={X1z    X2z;    X1z    X3z;    X2z    X3z};
        DLBL_test={X1_str X2_str; X1_str X3_str; X2_str X3_str};
    
    
      % ----- LOOP OVER TEST PREDICTORS
    
        for tt = 1 : Nlvl2
    
            it=it+1;
    
            % ---- Combine & compute metadata
            % Combine X & Y
              ALLD=cell(1,iy);                     ALLD(1:useN_b)=useD_base;
              ALLD(useN_b+1:ixd)=useD_test(tt,:);  ALLD{iy}=Yz_m; 
              PLTLBL_ALL{it}=[strjoin(DLBL_test(tt,:),'+') Dcdastr];
              disp([' --> Running test for: ' PLTLBL_ALL{it}])
    
                
            % ========== Global
    
              pval_g=zeros(1,Nfc);
              deltaR2_g=zeros(1,Nfc);
              for ii = 1 : Nfc
                  disp([' Global: ' ylbl{ii}])
                % Extract all X data with the FC data for this round
                  D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
                % Remove 0-valued edges
                  Dt=rmzeros(D); X=Dt(:,1:ixd); y=Dt(:,iy);
                % Modeling: Hiearchical
                  X_base=X(:,1:useN_b);
                  [deltaR2_g(ii),pval_g(ii)]=compare_models(y,X_base,X);
              end
            % Store outputs
              deltaR2_g_ALL(it,:)=deltaR2_g; 
              pval_g_ALL(it,:)=pval_g;
    
        
            % ========== Pairwise networks  
    
              pval_ntwk=zeros(Nntwk,Nntwk,Nfc);
              deltaR2_ntwk=zeros(Nntwk,Nntwk,Nfc);
              for ii = 1 : Nfc
                  disp([' Network: ' ylbl{ii}])
                % Extract all X data with the FC data for this round
                  D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
                  for iio = 1 : Nntwk
                      nio=pinfo.cis{iio};
                      for iii = iio : Nntwk
                          nii=pinfo.cis{iii};
                          ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(nio,nii); end 
                          Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
                        % Modeling: Hiearchical
                          X_base=X(:,1:useN_b);
                          [deltaR2_ntwk(iio,iii,ii),pval_ntwk(iio,iii,ii)]=compare_models(y,X_base,X);
                      end
                  end
              end
            % Store outputs
              deltaR2_ntwk_ALL(it,:,:,:)=deltaR2_ntwk; 
              pval_ntwk_ALL(it,:,:,:)=pval_ntwk;
    
    
           % ==========        Node-level    
    
             pval_node=zeros(Nnode,Nfc);
             deltaR2_node=zeros(Nnode,Nfc);
             for ii = 1 : Nfc
                  disp([' Node: ' ylbl{ii}])
                % Extract all X data with the FC data for this round
                  D=cell(1,iy); D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
                  for nn = 1 : Nnode
                      % Select subset of data (individual nodes)
                        ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(:,nn); end
                      % rm zeros
                        Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
                      % Modeling: Hiearchical
                        X_base=X(:,1:useN_b);
                        [deltaR2_node(nn,ii),pval_node(nn,ii)]=compare_models(y,X_base,X);
                  end
             end
            % Store outputs
              deltaR2_node_ALL(it,:,:)=deltaR2_node; 
              pval_node_ALL(it,:,:)=pval_node;
        end
    
      %% ======== LEVEL 3: all edge weights
      
        ixd=useN_b + 3;  iy=ixd+1;
        useD_test={X1z    X2z    X3z};
        DLBL_test={X1_str X2_str X3_str};
    
    
      % ----- LOOP OVER TEST PREDICTORS
    
        for tt = 1 : Nlvl3
    
            it=it+1;
    
            % ---- Combine & compute metadata
            % Combine X & Y
              ALLD=cell(1,iy);                     ALLD(1:useN_b)=useD_base;
              ALLD(useN_b+1:ixd)=useD_test(tt,:);  ALLD{iy}=Yz_m; 
              PLTLBL_ALL{it}=[strjoin(DLBL_test(tt,:),'+') Dcdastr];
              disp([' --> Running test for: ' PLTLBL_ALL{it}])
    
    
            % ========== Global  
                
              pval_g=zeros(1,Nfc);
              deltaR2_g=zeros(1,Nfc);
              for ii = 1 : Nfc
                  disp([' Global: ' ylbl{ii}])
                % Extract all X data with the FC data for this round
                  D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
                % Remove 0-valued edges
                  Dt=rmzeros(D); X=Dt(:,1:ixd); y=Dt(:,iy);
                % Modeling: Hiearchical
                  X_base=X(:,1:useN_b);
                  [deltaR2_g(ii),pval_g(ii)]=compare_models(y,X_base,X);
              end
            % Store outputs
              deltaR2_g_ALL(it,:)=deltaR2_g; 
              pval_g_ALL(it,:)=pval_g;
    
        
            % ========== Pairwise networks  
              
              pval_ntwk=zeros(Nntwk,Nntwk,Nfc);
              deltaR2_ntwk=zeros(Nntwk,Nntwk,Nfc);
              for ii = 1 : Nfc
                  disp([' Network: ' ylbl{ii}])
                % Extract all X data with the FC data for this round
                  D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
                  for iio = 1 : Nntwk
                      nio=pinfo.cis{iio};
                      for iii = iio : Nntwk
                          nii=pinfo.cis{iii};
                          ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(nio,nii); end 
                          Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
                        % Modeling: Hiearchical
                          X_base=X(:,1:useN_b);
                          [deltaR2_ntwk(iio,iii,ii),pval_ntwk(iio,iii,ii)]=compare_models(y,X_base,X);
                      end
                  end
              end
            % Store outputs
              deltaR2_ntwk_ALL(it,:,:,:)=deltaR2_ntwk; 
              pval_ntwk_ALL(it,:,:,:)=pval_ntwk;
    
    
            % ==========        Node-level  
    
              pval_node=zeros(Nnode,Nfc);
              deltaR2_node=zeros(Nnode,Nfc);
              for ii = 1 : Nfc
                  disp([' Node: ' ylbl{ii}])
                % Extract all X data with the FC data for this round
                  D=cell(1,iy); D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
                  for nn = 1 : Nnode
                      % Select subset of data (individual nodes)
                        ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(:,nn); end
                      % rm zeros
                        Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
                      % Modeling: Hiearchical
                        X_base=X(:,1:useN_b);
                        [deltaR2_node(nn,ii),pval_node(nn,ii)]=compare_models(y,X_base,X);
                  end
              end
            % Store outputs
              deltaR2_node_ALL(it,:,:)=deltaR2_node; 
              pval_node_ALL(it,:,:)=pval_node;
    
        end
    
          Dcdastr_ALL{ss}=Dcdastr;
    end
    
    
    
    %% ===== Visualization
    
      PLTLBL_ALL_simple=cellfun(@(c)  extractBefore(c, ' '),PLTLBL_ALL,'UniformOutput',0);
    
    
    %% ======= GLOBAL: 3-panel ΔR² plots (Levels 1–3)
    % tiledlayout BARplots
    
    str2 = ['Global, ' str]; 
    ppostmp1 = [-2559 237 670 560];  % thin
    % ppostmp = [-2559 237 586 560]; % too thin
    myfig([useCMLBL ': ' str2], ppostmp1);
    YLIM = ceil(max(deltaR2_g_ALL(:))*100)/100;
    
    t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    
    % ---------------- Level 1 ----------------
      cmap_tmp = cmap_pred(1:Nlvl1,:);
      D1 = deltaR2_g_ALL(1:Nlvl1,:); 
      D2 = pval_g_ALL(1:Nlvl1,:);
      D1(D2>0.05) = 0; 
      ax1=nexttile; B=bar(D1',0.8); font(fntsz,usefont);
      ylabel('\DeltaR^2','FontWeight','bold');
      for xx=1:Nlvl1; B(xx).FaceColor=cmap_tmp(xx,:); end
      grid on; set(gca,'XGrid','off','TickLength',[0 0]); ylim([0 YLIM]);
    
    % ---------------- Level 2 ----------------
      cmap_tmp = gray(Nlvl2+2);
      D1 = deltaR2_g_ALL(Nlvl1+1:Nlvl1+Nlvl2,:); 
      D2 = pval_g_ALL(Nlvl1+1:Nlvl1+Nlvl2,:);
      D1(D2>0.05) = 0; 
      ax2=nexttile; B=bar(D1',0.8); font(fntsz,usefont);
      ylabel('\DeltaR^2','FontWeight','bold');
      for xx=1:Nlvl2; B(xx).FaceColor=cmap_tmp(xx+1,:); end
      grid on; set(gca,'XGrid','off','TickLength',[0 0]); ylim([0 YLIM]);
    
    % ---------------- Level 3 ----------------
      cmap_tmp = [0 0 0];
      D1 = deltaR2_g_ALL(Ntot,:); 
      D2 = pval_g_ALL(Ntot,:);
      D1(D2>0.05) = 0;  
      ax3=nexttile; B=bar(D1',0.5); font(fntsz,usefont);
      ylabel('\DeltaR^2','FontWeight','bold');
      B.FaceColor=cmap_tmp;
      grid on; set(gca,'XGrid','off'); ylim([0 YLIM]);
      set(gca,'XTickLabels',ylbl,'XTickLabelRotation',35,'TickLength',[0 0]);
    
      set([ax1,ax2,ax3],'XLim',[0.5 Nfc+0.5],'XTick',1:Nfc,'XTickLabel',ylbl);
      set([ax1,ax2],'XTickLabel','');
    
    % Save
      svstr='1_global';
      set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '_' useCMLBL '.pdf'],'ContentType','vector');
    
    
    %% ======= Network: 3-panel ΔR² plots (Levels 1–3)
    % tiledlayout BOXplots
    
    % Preallocation
      Delta_ntwk=cell(Ntot,Nfc); Pval_ntwk=cell(Ntot,Nfc);                                         % 7 levels x Nfc y's
      for lvl = 1:Ntot
            for yy = 1:Nfc
                Dmat=squeeze(deltaR2_ntwk_ALL(lvl,:,:,yy)); Pmat=squeeze(pval_ntwk_ALL(lvl,:,:,yy)); % extract data
                mask=triu(true(size(Dmat))); vals=Dmat(mask); pvals=Pmat(mask);                      % take only upper tri + main diag
                Delta_ntwk{lvl,yy}=vals(:); Pval_ntwk{lvl,yy}=pvals(:);                              % store
            end
      end
    
    % ===== PLOT
      str2=['Network ΔR², ' str]; ppostmp2 = [-2559 237 832 560];  % thick
      myfig([useCMLBL ': ' str2],ppostmp2);
      YLIM=ceil(max(deltaR2_ntwk_ALL(:))*100)/100;
      boxopt={'BoxWidth',0.2,'MarkerStyle','*','MarkerSize',8,'BoxFaceAlpha',.9,'BoxMedianLineColor','k','LineWidth',1.5};
      boxopt2={'BoxWidth',0.5,'MarkerStyle','*','MarkerSize',8,'BoxFaceAlpha',.9,'BoxMedianLineColor','w','BoxFaceColor','k','MarkerColor','k','LineWidth',1.5};
      w=0.5; offsets=linspace(-w/2,w/2,Nlvl1);                                     % evenly spaced offsets for levels 1 & 2
      t=tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    
    % ---- Level 1 ----
    ax1 = nexttile; hold on;
    for pr = 1:Nlvl1
        for yy = 1:Nfc
            data  = Delta_ntwk{pr,yy};
            pvals = Pval_ntwk{pr,yy};
            data(pvals > 0.05) = [];
    
            x = yy + offsets(pr);      % apply offset for each predictor
            boxchart(x*ones(size(data)),data,'BoxFaceColor',cmap_pred(pr,:),'MarkerColor',cmap_pred(pr,:),boxopt{:});
        end
    end
    ylabel('\DeltaR^2','FontWeight','bold'); font(fntsz,usefont);
    ylim([0 YLIM]); grid on; set(gca,'XGrid','off','TickLength',[0 0]);
    
    % ---- Level 2 ----
    cmap_tmp = gray(Nlvl2+2);
    ax2 = nexttile; hold on;
    for pr = 1:Nlvl2
        for yy = 1:Nfc
            data  = Delta_ntwk{Nlvl1+pr,yy};
            pvals = Pval_ntwk{Nlvl1+pr,yy};
            data(pvals > 0.05) = [];
    
            x = yy + offsets(pr);
            boxchart(x*ones(size(data)),data,boxopt{:},'BoxFaceColor',cmap_tmp(pr+1,:),'MarkerColor',cmap_tmp(pr+1,:));
        end
    end
    ylabel('\DeltaR^2','FontWeight','bold'); font(fntsz,usefont);
    ylim([0 YLIM]); grid on; set(gca,'XGrid','off','TickLength',[0 0]);
    
    % ---- Level 3 ----
    ax3 = nexttile; hold on;
    for yy = 1:Nfc
        data = Delta_ntwk{Ntot,yy};
        pvals = Pval_ntwk{Ntot,yy};
        data(pvals > 0.05) = [];
    
        boxchart(yy*ones(size(data)),data,boxopt2{:});
    end
    ylabel('\DeltaR^2','FontWeight','bold'); font(fntsz,usefont);
    ylim([0 YLIM]); grid on; set(gca,'XGrid','off','TickLength',[0 0]);
    
    % ---- Axis formatting
      set([ax1,ax2,ax3],'XLim',[0.5 Nfc+0.5],'XTick',1:Nfc,'Box','on');
      set([ax1,ax2],'XTickLabel','');
      set(ax3,'XTickLabel',ylbl,'XTickLabelRotation',35);
    
    
    % Save
      svstr='2_ntwk';
      set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '_' useCMLBL '.pdf'],'ContentType','vector');
    
    
    
    % % % %% =======  Network-Level: Matrix Plots of ΔR²
    % % %   axopt={'XTick',1:Nntwk,'XTickLabel',LBL_ntwk,'XTickLabelRotation',40, ...
    % % %          'YTick',1:Nntwk,'YTickLabel',LBL_ntwk,'YTickLabelRotation',40,'FontSize',fntsz};
    % % %   meshopt={'EdgeColor','k','LineWidth',1.5};
    % % %   ppostmp2=getPlotPos([-2559 534 315 263],Ntot,0);
    % % % % --- ΔR²
    % % %   ii=1;       % switch for FC
    % % %   xx=[1 2 3 7]; % switch for tests
    % % %   D=squeeze(deltaR2_ntwk_ALL(xx,:,:,ii)); D2=squeeze(pval_ntwk_ALL(xx,:,:,ii));
    % % %   D(D<0)=0; D(D2>0.05)=0; TMPSTR='ΔR²: '; usecmap=cm3; 
    % % %   TMPSTR2=[TMPSTR useCMLBL ' ' ylbl{ii}];  TMPLBLS=PLTLBL_ALL_simple(xx);
    % % %   CLIM=[0 ceil(prctile(D(:),97)*100)/100]; %CLIM=[0 ceil(max(D(:))*10)/10];
    % % %   for xxx=1:length(xx)
    % % %     % Matrix with blank upper tri
    % % %       D2=squeeze(D(xxx,:,:)); tmask=tril(true(size(D2))); myfig([TMPSTR2 ', ' TMPLBLS{xxx}],ppostmp2(xxx,:));
    % % %       imagesc(D2','AlphaData',tmask); colormap(usecmap); colorbar; axis square; clim(CLIM); font(fntsz,usefont);
    % % %       set(gca,axopt{:}); font(fntsz,usefont); hold on; title(TMPLBLS{xxx}); 
    % % %       for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
    % % %       set(gca,'XTickLabel','','YTickLabel',''); title ''; axis off;
    % % %       % colorbar off;
    % % %   end
    
    
    
    %% =======  Node-Level: Surface Plots of ΔR²
    
    % ----- ΔR²
      xx=[1 2 3 7]; % switch for tests 
    
    
      for ii=[1 5 8]     % switch for FC 
          D=squeeze(deltaR2_node_ALL(xx,:,ii)); D2=squeeze(pval_node_ALL(xx,:,ii)); TMPLBLS=PLTLBL_ALL_simple(xx);
          D(D<0)=0; D(D2>0.05)=0; TMPSTR='ΔR²: '; usecmap=[.75 .75 .75; cm4]; 
          TMPSTR2=[TMPSTR useCMLBL ' ' ylbl{ii}]; CLIM=[0 ceil(prctile(D(:),97)*100)/100]; %CLIM=[0 ceil(max(D(:))*10)/10];
        
        % % % % Surface Plot (all together)
        % % %   Hcx=plot_conn_surf(D',pinfo,'cortex',TMPLBLS,TMPSTR2); setsurf(Hcx,CLIM,usecmap);
        % % % % Save
        % % %   svstr=['3_node_' num2str(CC) '-' useCMLBL '_' num2str(ii) '-' strrep(ylbl{ii},'_','')];
        % % %   set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '.pdf'],'ContentType','vector');
        
        % Surface Plot (individual)
        for dd=1:4
          Hcx=plot_conn_surf(D(dd,:)',pinfo,'cortex',TMPLBLS(dd),TMPSTR2); setsurf(Hcx,CLIM,usecmap);
          % Save
          svstr=['3_node_' num2str(CC) '-' useCMLBL '_' num2str(ii) '-' strrep(ylbl{ii},'_','') '_' num2str(dd) '-' strrep(TMPLBLS{dd},'+','-')];
          set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '.pdf'],'ContentType','vector');
        end
      end

end



% -- Corr with SAaxis
  dx=SA_map; dxlbl='S-A map'; savecorr=zeros(1,length(xx)); ppostmp=[-2544 394 318 403];
  for jj=1:length(xx)
      dy=D(jj,:); dy=dy'; savecorr(jj)=corr(dx(dy~=0),dy(dy~=0));
  end
% Bar plot
  myfig([TMPSTR2 ', vs SAaxis'],ppostmp); B=bar(TMPLBLS,savecorr,0.6); % [-1115 359 467 330]
  ylabel('Pearson''s r'); font(fntsz,usefont); % title([ 'R^2 vs ' dxlbl]);
  grid on; set(gca,'XGrid','off','XMinorGrid','off');
  B.CData=CLRS.carolinablue; B.FaceColor='flat';

  B.CData=[cmap_pred(1:3,:); 0 0 0]; ylim([-.59 .59]); 
  
  % set(gca,'XTickLabels',''); set(gcf,'Position',[-2559 547 199 250]); ylabel '';



% -- Corr with SAaxis: loop over all FC
  xx=[1 2 3 7]; % switch for tests 
  ppostmp=getPlotPos([-2559 547 199 250],Nfc,0);

  for ii=1:Nfc     % switch for FC 
      D=squeeze(deltaR2_node_ALL(xx,:,ii)); D2=squeeze(pval_node_ALL(xx,:,ii)); TMPLBLS=PLTLBL_ALL_simple(xx);
      D(D<0)=0; D(D2>0.05)=0; TMPSTR='ΔR²: ';
      TMPSTR2=[TMPSTR useCMLBL ' ' ylbl{ii}]; CLIM=[0 ceil(max(D(:))*10)/10];
    
      dx=SA_map; dxlbl='S-A map'; savecorr=zeros(1,length(xx)); 
      for jj=1:length(xx)
          dy=D(jj,:); dy=dy'; savecorr(jj)=corr(dx(dy~=0),dy(dy~=0));
      end
    % Bar plot
      myfig([TMPSTR2 ', vs SAaxis'],ppostmp(ii,:)); B=bar(TMPLBLS,savecorr,0.6);
      ylabel('Pearson''s r'); font(fntsz,usefont); 
      grid on; set(gca,'XGrid','off','XMinorGrid','off');
      B.CData=CLRS.carolinablue; B.FaceColor='flat';
    
      B.CData=[cmap_pred(1:3,:); 0 0 0]; ylim([-.59 .59]); 
      
      set(gca,'XTickLabels',''); ylabel '';
  end

% -------------------------------------------------------------------------




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% =======   Modeling 1b: Blockwise RSN FC prediction
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% All tests run within a single communication model


% % % % % % % Inputs:
% % % % % % % rsn_labels (N×1), rsn_names (1×R), rsn_tier (R×1 cellstr/categorical: 'unimodal'/'transmodal')
% % % % % % 
% % % % % % 
% % % % % % % Build FC_struct
% % % % % %   FC_struct = make_fc_struct(Y_m, ylbl);
% % % % % % 
% % % % % % % Build Comm_struct
% % % % % % useCM=[ALLCM {ALLCM_delay}]; 
% % % % % % useCMLBL=[CMLBL {CMLBL_delay}]; 
% % % % % % dsLabels = {'caliber','MTsat','g-ratio','delay'};
% % % % % %   Comm_struct = make_comm_struct(useCM, useCMLBL, dsLabels, CMLBL_simple);
% % % % % % 
% % % % % % 
% % % % % % opts = struct('nBoot', 500, 'standardize_predictors', true, 'verbose', true);
% % % % % % out = blockwise_rsn_predict_fc(FC_struct, Comm_struct, rsn_labels, rsn_names, rsn_tier, opts);
% % % % % % 
% % % % % % % Inspect results:
% % % % % % out.results(:,{'predictor','family','weight_type','dR2','beta_zxUW','p_zxUW','beta_zxTB','p_zxTB'})
% % % % % % out.replication
% % % 
% % % CC=5;    % switch for communication models
% % % 
% % % 
% % % % --- Prep inputs ---
% % % EDz=nzzscore(EDnrm,0);
% % % base = struct('Xb', ALLCM_bin{CC}, ...
% % %               'Xc', ALLCM{1}{CC}, ...
% % %               'Xd', EDz, ...
% % %               'labels_base', {'binary','caliber','ED'});
% % % 
% % % tests = struct('X1', ALLCM{2}{CC}, ...
% % %                'X2', ALLCM{3}{CC}, ...
% % %                'X3', ALLCM_delay{CC}, ...
% % %                'labels_tests', { {'MTsat','g-ratio','delay'} });
% % % 
% % % optsA = struct('alpha',0.05,'use_lower_only',true,'label_model',CMLBL_simple{CC},'include_pairs',false);
% % % 
% % % pinfo.clabels_short=LBL_ntwk;
% % % 
% % % % --- Stage A ---
% % % [Delta_block, P_block, stats] = fit_BICS_stageA(Yz_m, base, tests, pinfo, optsA);
% % % stats.model_labels{end}='all';
% % % 
% % % % --- Stage B ---
% % % optsB = struct('nperm',200,'eps_scale',0.05,'use_lower_only',true,'rng_seed',13, ...
% % %                'block_map',stats.block_map,'perm_mode','Xjoint','includePairs',stats.include_pairs); %'reg_mode','fitlm'
% % % [SI_block, elasticity, pert] = fit_BICS_stageB(Yz_m, base, tests, pinfo, Delta_block, optsB);
% % % 
% % % % --- Viz ---
% % % cmaptmp=linspecer(Nfc);
% % % viz_opts = struct('fc_agg','mean','cmap_delta',cm,'cmap_SI',cm2,    ...
% % %                   'cmap_fc', cmaptmp,'figpos',[-2559 384 1400 413], ...
% % %                   'elasticity',elasticity,'elasticity_norm','by_delta', ...
% % %                   'elasticity_norm_floor',1e-6,'elasticity_as_percent',true);
% % % viz_BICS(Delta_block, SI_block, stats, ylbl, pinfo, viz_opts);






% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% =======   Modeling 1c: DUAL Blockwise RSN FC prediction
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% All tests designed to include 1 routing & 1 diffusion model together
% Interactions added as well

% --- Setup
  CC                    = [1 5];
  CMstr                 = strjoin(CMLBL_simple(CC),'-');
  ALLY                  = Yz_m;
  EDz                   = nzzscore(EDnrm,0);
  ALLX                  = [ALLCM {ALLCM_delay} {ALLCM_bin} {EDz}]; 
  pinfo.clabels_short   = LBL_ntwk;
  uselower              = 1; 
  elast_mode            = 'zrow';
  FCAGG                 = 'median';

% DEBUGGING
  % % % optsA.debug = struct('probe',true,'probe_blocks',false,'probe_nodes',false,'trace_menus',true,'save','dbg_stageA.mat');
  % % % optsB.debug = struct('probe_B',true,'save','dbg_stageB.mat');
  % % % tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
  % % %              'label_comm_models',strjoin(CMLBL_simple(CC),'-'), ...
  % % %              'modes', 'all', 'levels', [1 2],                   ...
  % % %              'interaction', 'caliber', 'interact_at', 'both');
  % % % 
  % % %              % 'interaction', 'none');
  % % % [Delta_all,~,~] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, optsA);
  % % % [~,~,~] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, optsB);

% Run Opts
  opts = struct( ...
            'stageA', struct( ...
                        'use_lower_only', uselower, ...
                        'alpha', 0.05 ...                                   % intended for FDR (add later)
                            ), ...
            'stageB', struct( ...
                        'nperm', 200, ...
                        'eps_scale', 0.05, ...
                        'use_lower_only', uselower, ...
                        'rng_seed', 13, ...
                        'perm_mode', 'yresid', ...
                        'elasticity_mode', elast_mode ...                   % 'xscale' | 'zadd'
                        ), ...
            'viz', struct( ...                                              % main results
                        'figpos1', [-2559 54 563 743], ...                  % heatmaps (M=1 case)
                        'labelON', true, ...                                % turn some plot labels on/off
                        'fc_agg', FCAGG, ...                                % 'mean'|'median'
                        'elasticity', [], ...                               % for Figure 3 style panels
                        'elasticity_norm', 'none', ...                      % 'none'|'by_delta'
                        'elasticity_norm_floor', 1e-6, ...
                        'elasticity_as_percent', true, ...
                        'ribbons_xdim', 'bands' ...                         % 'models' | 'bands'
                        ), ...
            'vizE', struct( ...                                             % elastcitiy results
                        'figpos1', [-2559 54 321 743], ...                  % heatmaps (M=1 case)
                        'fc_agg', FCAGG, ...                                % 'mean'|'median'
                        'labelON', true, ...                                % turn some plot labels on/off
                        'elasticity_mode', elast_mode, ...                  % 'xscale' | 'zadd'
                        'norm', 'none', ...                                 % 'none'|'by_delta'
                        'norm_floor', 1e-6, ...
                        'ribbons_xdim', 'bands' ...                         % 'models' | 'bands'
                        ), ...
            'vizN', struct( ...                                             % node results
                        'figpos_bars', [-936 552 936 245], ...              % M=1 case
                        'corr_type', 'Spearman', ...                        % 'mean'|'median'
                        'absval', false, ...                                % 'xscale' | 'zadd'
                        'm_eq1_colors', [], ...                             % 'none'|'by_delta'
                        'bands', 1:Nfc, ...                                 % numeric
                        'labelON', true, ...                                % turn some plot labels on/off
                        'make_surfaces', true, ...
                        'surf_bands', [], ...
                        'surf_level', 'L2_both',...
                        'surf_title', 'deltaR2', ...
                        'surf_cmax', [], ...                                % auto (95th pct of |ΔR²|)
                        'surf_cmap', [.7 .7 .7; cm4] ... 
                        ) ...
                );

% -------------------------------------------------------------------------
% =====  Base (no interactions) =====
% -------------------------------------------------------------------------

tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
               'label_comm_models',CMstr, ...
               'modes', 'all', 'levels', [1 2],                   ...
               'interaction', 'none', 'effect_mode', 'incremental');


% STAGE A: Compute ΔR² per RSN block, model spec, and FC target.
  [Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);

% STAGE B: Sensitivity and elasticity testing
  [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts.stageB);

% Visualization: deltaR2 & SI
  figs = bicsdual.viz_BICS_dual(Delta_all, SI_all, stats_all, ylbl, pinfo, opts.viz);
  set(figs.ribbons_Meq1,'Position',[-2559 505 729 292]); AX=findall(figs.ribbons_Meq1,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);            % line
  set(figs.heatmaps_Meq1,'Position',[-2559 54  563 743]);                   % heatmaps
  
% Visualization: elasticity
  figsE = bicsdual.viz_BICS_dual_elast(elasticity_all, stats_all, ylbl, pinfo, opts.vizE);
  set(figsE.elast_ribbons_Meq1,'Position',[-728 505 729 292]); title ''; xlabel ''; xtickangle(30);  % line
  set(figsE.elast_heatmaps_Meq1,'Position',[-320 54 321 743])                                        % heatmaps
  

% Visualization: surf & S-A axis corr
  figsN = bicsdual.sa_corr_nodal(stats_all, SA_map, ylbl, pinfo, opts.vizN);
  set(figsN.fig_surf(1),'Position',[-1.7778   -0.0156    0.9000    0.9000])
  set(figsN.fig_surf(2),'Position',[-1.1444   -0.0156    0.9000    0.9000])
  set(figsN.fig_bars,'Position',[-2559 -539 729 245]); AX=findall(figsN.fig_bars,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);  % line


% -------------------------------------------------------------------------
% ==========          CALIBER interactions  (L1 and L2)          ==========
% -------------------------------------------------------------------------

tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
               'label_comm_models',CMstr, ...
               'modes', 'all', 'levels', [1 2],                   ...
               'interaction', 'caliber', 'interact_at', 'both', 'effect_mode', 'incremental');

% STAGE A: Compute ΔR² per RSN block, model spec, and FC target.
  [Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);

% STAGE B: Sensitivity and elasticity testing
  [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts.stageB);

% Visualization: deltaR2 & SI
  figs = bicsdual.viz_BICS_dual(Delta_all, SI_all, stats_all, ylbl, pinfo, opts.viz);
  set(figs.ribbons_Meq1,'Position',[-2559 505 729 292]); AX=findall(figs.ribbons_Meq1,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);                               % line  
  set(figs.heatmaps_Meq1,'Position',[-2559 54  563 743]);                                      % heatmaps
  
% Visualization: elasticity
  figsE = bicsdual.viz_BICS_dual_elast(elasticity_all, stats_all, ylbl, pinfo, opts.vizE);
  set(figsE.elast_ribbons_Meq1,'Position',[-728 505 729 292]); title ''; xlabel ''; xtickangle(30);  % line
  set(figsE.elast_heatmaps_Meq1,'Position',[-320 54 321 743])                                        % heatmaps
  

% Visualization: surf & S-A axis corr
  figsN = bicsdual.sa_corr_nodal(stats_all, SA_map, ylbl, pinfo, opts.vizN);
  set(figsN.fig_surf(1),'Position',[-1.7778   -0.0156    0.9000    0.9000])
  set(figsN.fig_surf(2),'Position',[-1.1444   -0.0156    0.9000    0.9000])
  set(figsN.fig_bars,'Position',[-2559 -539 729 245]); AX=findall(figsN.fig_bars,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);  % line


% -------------------------------------------------------------------------
% ==========             ED interactions  (L1 and L2)            ==========
% -------------------------------------------------------------------------

tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
               'label_comm_models',CMstr, ...
               'modes', 'all', 'levels', [1 2],                   ...
               'interaction', 'ed', 'interact_at', 'both', 'effect_mode', 'incremental');

% STAGE A: Compute ΔR² per RSN block, model spec, and FC target.
  [Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);

% STAGE B: Sensitivity and elasticity testing
  [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts.stageB);

% Visualization: deltaR2 & SI
  figs = bicsdual.viz_BICS_dual(Delta_all, SI_all, stats_all, ylbl, pinfo, opts.viz);
  set(figs.ribbons_Meq1,'Position',[-2559 505 729 292]); AX=findall(figs.ribbons_Meq1,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);                               % line  
  set(figs.heatmaps_Meq1,'Position',[-2559 54  563 743]);                                      % heatmaps
  
% Visualization: elasticity
  figsE = bicsdual.viz_BICS_dual_elast(elasticity_all, stats_all, ylbl, pinfo, opts.vizE);
  set(figsE.elast_ribbons_Meq1,'Position',[-728 505 729 292]); title ''; xlabel ''; xtickangle(30);  % line
  set(figsE.elast_heatmaps_Meq1,'Position',[-320 54 321 743])                                        % heatmaps
  

% Visualization: surf & S-A axis corr
  figsN = bicsdual.sa_corr_nodal(stats_all, SA_map, ylbl, pinfo, opts.vizN);
  set(figsN.fig_surf(1),'Position',[-1.7778   -0.0156    0.9000    0.9000])
  set(figsN.fig_surf(2),'Position',[-1.1444   -0.0156    0.9000    0.9000])
  set(figsN.fig_bars,'Position',[-2559 -539 729 245]); AX=findall(figsN.fig_bars,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);  % line


% -------------------------------------------------------------------------
% ==========            BOTH interactions  (L1 and L2)           ==========
% -------------------------------------------------------------------------

tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
               'label_comm_models',CMstr, ...
               'modes', 'all', 'levels', [1 2],                   ...
               'interaction', 'both', 'interact_at', 'both', 'effect_mode', 'incremental');

% STAGE A: Compute ΔR² per RSN block, model spec, and FC target.
  [Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);

% STAGE B: Sensitivity and elasticity testing
  [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts.stageB);

% Visualization: deltaR2 & SI
  figs = bicsdual.viz_BICS_dual(Delta_all, SI_all, stats_all, ylbl, pinfo, opts.viz);
  set(figs.ribbons_Meq1,'Position',[-2559 505 729 292]); AX=findall(figs.ribbons_Meq1,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);                               % line  
  set(figs.heatmaps_Meq1,'Position',[-2559 54  563 743]);                                      % heatmaps
  
% Visualization: elasticity
  figsE = bicsdual.viz_BICS_dual_elast(elasticity_all, stats_all, ylbl, pinfo, opts.vizE);
  set(figsE.elast_ribbons_Meq1,'Position',[-728 505 729 292]); title ''; xlabel ''; xtickangle(30);  % line
  set(figsE.elast_heatmaps_Meq1,'Position',[-320 54 321 743])                                        % heatmaps
  

% Visualization: surf & S-A axis corr
  figsN = bicsdual.sa_corr_nodal(stats_all, SA_map, ylbl, pinfo, opts.vizN);
  set(figsN.fig_surf(1),'Position',[-1.7778   -0.0156    0.9000    0.9000])
  set(figsN.fig_surf(2),'Position',[-1.1444   -0.0156    0.9000    0.9000])
  set(figsN.fig_bars,'Position',[-2559 -539 729 245]); AX=findall(figsN.fig_bars,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);  % line



% -------------------------------------------------------------------------
% =======     Total contribution (main effects + interactions)      =======
% -------------------------------------------------------------------------

tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
               'label_comm_models',CMstr, ...
               'modes', 'all', 'levels', [1 2],                   ...
               'interaction', 'both', 'interact_at', 'both', 'effect_mode', 'total');

% STAGE A: Compute ΔR² per RSN block, model spec, and FC target.
  [Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);

% STAGE B: Sensitivity and elasticity testing
  [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts.stageB);

% Visualization: deltaR2 & SI
  figs = bicsdual.viz_BICS_dual(Delta_all, SI_all, stats_all, ylbl, pinfo, opts.viz);
  set(figs.ribbons_Meq1,'Position',[-2559 505 729 292]); AX=findall(figs.ribbons_Meq1,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);                               % line  
  set(figs.heatmaps_Meq1,'Position',[-2559 54  563 743]);                                      % heatmaps
  
% Visualization: elasticity
  figsE = bicsdual.viz_BICS_dual_elast(elasticity_all, stats_all, ylbl, pinfo, opts.vizE);
  set(figsE.elast_ribbons_Meq1,'Position',[-728 505 729 292]); title ''; xlabel ''; xtickangle(30);  % line
  set(figsE.elast_heatmaps_Meq1,'Position',[-320 54 321 743])                                        % heatmaps
  

% Visualization: surf & S-A axis corr
  figsN = bicsdual.sa_corr_nodal(stats_all, SA_map, ylbl, pinfo, opts.vizN);
  set(figsN.fig_surf(1),'Position',[-1.7778   -0.0156    0.9000    0.9000])
  set(figsN.fig_surf(2),'Position',[-1.1444   -0.0156    0.9000    0.9000])
  set(figsN.fig_bars,'Position',[-2559 -539 729 245]); AX=findall(figsN.fig_bars,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);  % line



% -------------------------------------------------------------------------
%% ---------------       Load & test saved results     --------------------
% -------------------------------------------------------------------------


% --- Setup
  ALLY                  = Yz_m;
  EDz                   = nzzscore(EDnrm,0);
  ALLX                  = [ALLCM {ALLCM_delay} {ALLCM_bin} {EDz}]; 
  pinfo.clabels_short   = LBL_ntwk;


for xx=1:4;
close all;
% ---- Overall Setup
  testset               = {[1 5] [1 6] [2 5] [2 6]};                        % SPE-CMY  SPE-DE  NE-CMY  NE-DE
  CC                    = testset{xx};
  CMstr                 = strjoin(CMLBL_simple(CC),'-');

% Load data
  svdir=[locDeriv '/x_dual_BICS/grp_' parc];   
  % load([svdir '/results_' CMstr '.mat'])
  load([svdir '/results_total_tmp' CMstr '.mat'])   % version that includes total contributions


% Extract metadata
  opts              = S.meta.opts;
  label_comm_models = S.meta.label_comm_models;
  myelin_preds      = S.meta.myelin_predictors{:};

% ===== Do stuff with a single 

% -- Switch between interactions
  ii                = 3;
  label_interact    = {'int_none' 'int_caliber' 'int_ed' 'int_both' 'total'};
  intstr            = label_interact{ii}; 

  Delta_all         = S.(intstr).Delta_all;
  P_all             = S.(intstr).P_all;
  stats_all         = S.(intstr).stats_all;
  SI_all            = S.(intstr).SI_all;
  elasticity_all    = S.(intstr).elasticity_all;
  pert_all          = S.(intstr).pert_all;

  if ii<5
      stats_all.meta.effect_mode='incremental';
  else
      stats_all.meta.effect_mode='total';
  end

% % % %  --- Main visualizations
% % % % Visualization: deltaR2 & SI
% % %   figs = bicsdual.viz_BICS_dual(Delta_all, SI_all, stats_all, ylbl, pinfo, opts.viz);
% % %   set(figs.ribbons_Meq1,'Position',[-2559 505 729 292]); AX=findall(figs.ribbons_Meq1,'type','axes');
% % %   title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);                               % line  
% % %   set(figs.heatmaps_Meq1,'Position',[-2559 54  563 743]);                                      % heatmaps
% % %   FIG1=figs.heatmaps_Meq1;
% % %   % AX = findall(FIG4, 'type', 'axes');
% % % %   CLIM=[0 0.23]; 
% % % % AX=findall(FIG1,'type','axes'); AX(2).CLim=CLIM; AX(4).CLim=CLIM; AX(6).CLim=CLIM;
% % % % AX=findall(FIG2,'type','axes'); AX(2).CLim=CLIM; AX(4).CLim=CLIM; AX(6).CLim=CLIM;
% % % % AX=findall(FIG3,'type','axes'); AX(2).CLim=CLIM; AX(4).CLim=CLIM; AX(6).CLim=CLIM;
% % % % AX=findall(FIG4,'type','axes'); AX(2).CLim=CLIM; AX(4).CLim=CLIM; AX(6).CLim=CLIM;


% % % % ----- Visualization: deltaR2 & SI for all FC: match clims
% % %   FIG1 = cell(1,3); FIG2 = cell(1,3); AX1 = cell(1,3); AX2 = cell(1,3);
% % %   for k=1:3
% % %       opts.viz.usek=k;                                                          % Sets which level to use (e.g., L1-route, L2-both)
% % %       opts.viz.matchClim=true;                                                 % common clim across FC
% % %       figs = bicsdual.viz_BICS_dual_allFC(Delta_all, SI_all, stats_all, ylbl, pinfo, opts.viz);
% % %       % FIG1=figs.heatmaps1; FIG2=figs.heatmaps2;
% % %       FIG1{k}=figs.heatmaps1; FIG2{k}=figs.heatmaps2;
% % %       AX1{k} = findall(FIG1{k},'type','axes'); AX2{k} = findall(FIG2{k},'type','axes');
% % %   end
% % % % % % % 1 CLIM for all
% % % % % %   CLIM=[0 0.145]; 
% % % % % %   AX=findall(FIG1,'type','axes'); AX(2).CLim=CLIM; AX(4).CLim=CLIM; AX(6).CLim=CLIM; AX(8).CLim=CLIM;
% % % % % %   AX=findall(FIG2,'type','axes'); AX(2).CLim=CLIM; AX(4).CLim=CLIM; AX(6).CLim=CLIM; AX(8).CLim=CLIM;
% % % % CLIMS set by FC: FIG1
% % %   % CLIM_theta=[0 0.2]; CLIM_delta=[0 0.1]; CLIM_boldout=[0 0.2]; CLIM_boldin=[0 0.15];   % Main effects
% % %   % CLIM_theta=[0 0.05]; CLIM_delta=[0 0.04]; CLIM_boldout=[0 0.07]; CLIM_boldin=[0 0.11];   % int-caliber
% % %   % CLIM_theta=[0 0.07]; CLIM_delta=[0 0.05]; CLIM_boldout=[0 0.07]; CLIM_boldin=[0 0.11];   % int-ED
% % %   % CLIM_theta=[0 0.25]; CLIM_delta=[0 0.15]; CLIM_boldout=[0 0.225]; CLIM_boldin=[0 0.225];   % total
% % %   CLIM_theta=[0 0.25]; CLIM_delta=[0 0.15]; CLIM_boldout=[0 0.225]; CLIM_boldin=[0 0.3];   % matched-clims
% % %   for k=1:3, AX1{k}(2).CLim=CLIM_theta; end
% % %   for k=1:3, AX1{k}(4).CLim=CLIM_delta; end
% % %   for k=1:3, AX1{k}(6).CLim=CLIM_boldout; end
% % %   for k=1:3, AX1{k}(8).CLim=CLIM_boldin; end
% % % % CLIMS set by FC: FIG2
% % %   % CLIM_gammahi=[0 0.25]; CLIM_gammalo=[0 0.2]; CLIM_beta=[0 0.2]; CLIM_alpha=[0 0.28];  % Main effects
% % %   % CLIM_gammahi=[0 0.1]; CLIM_gammalo=[0 0.1]; CLIM_beta=[0 0.08]; CLIM_alpha=[0 0.07];  % int-caliber
% % %   % CLIM_gammahi=[0 0.14]; CLIM_gammalo=[0 0.18]; CLIM_beta=[0 0.09]; CLIM_alpha=[0 0.11];  % int-ED
% % %   % CLIM_gammahi=[0 0.325]; CLIM_gammalo=[0 0.275]; CLIM_beta=[0 0.3]; CLIM_alpha=[0 0.35];  % total
% % %   CLIM_gammahi=[0 0.3]; CLIM_gammalo=[0 0.275]; CLIM_beta=[0 0.3]; CLIM_alpha=[0 0.3];  % matched-clims
% % %   for k=1:3, AX2{k}(2).CLim=CLIM_gammahi; end
% % %   for k=1:3, AX2{k}(4).CLim=CLIM_gammalo; end
% % %   for k=1:3, AX2{k}(6).CLim=CLIM_beta; end
% % %   for k=1:3, AX2{k}(8).CLim=CLIM_alpha; end
% % % 
% % %   svstr=num2str(xx);
% % %   set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '_' CMstr '_1_route_1.pdf'],'ContentType','vector');
% % %   set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '_' CMstr '_1_route_2.pdf'],'ContentType','vector');
% % %   set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '_' CMstr '_2_diff_1.pdf'],'ContentType','vector');
% % %   set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '_' CMstr '_2_diff_2.pdf'],'ContentType','vector');
% % %   set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '_' CMstr '_3_both_1.pdf'],'ContentType','vector');
% % %   set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr '_' CMstr '_3_both_2.pdf'],'ContentType','vector');



% ----- Visualization: deltaR2 & SI for all FC: heterogeneous clims
  FIG1 = cell(1,3); FIG2 = cell(1,3); AX1 = cell(1,3); AX2 = cell(1,3); 
  for k=1:3
      opts.viz.usek=k; opts.viz.matchClim=false;
      figs = bicsdual.viz_BICS_dual_allFC(Delta_all, SI_all, stats_all, ylbl, pinfo, opts.viz);
      FIG1{k}=figs.heatmaps1; FIG2{k}=figs.heatmaps2;
      AX1{k} = findall(FIG1{k},'type','axes'); AX2{k} = findall(FIG2{k},'type','axes');
  end
  for iax = 2:2:8
        D=cell2mat(cellfun(@(c) c(iax).CLim,AX1,'UniformOutput',0)); CLIM=prctile(D(D~=0),75); 
        for k=1:3, AX1{k}(iax).CLim=[0 CLIM]; end
        D=cell2mat(cellfun(@(c) c(iax).CLim,AX2,'UniformOutput',0)); CLIM=prctile(D(D~=0),75); 
        for k=1:3, AX2{k}(iax).CLim=[0 CLIM]; end
  end

  svstr=num2str(xx);
  set(gcf,'Renderer','painters'); exportgraphics(FIG1{1},['~/Downloads/' svstr '_' CMstr '_1_route_FC1.pdf'],'ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(FIG1{2},['~/Downloads/' svstr '_' CMstr '_2_diff_FC1.pdf'],'ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(FIG1{3},['~/Downloads/' svstr '_' CMstr '_3_both_FC1.pdf'],'ContentType','vector');

  set(gcf,'Renderer','painters'); exportgraphics(FIG2{1},['~/Downloads/' svstr '_' CMstr '_1_route_FC2.pdf'],'ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(FIG2{2},['~/Downloads/' svstr '_' CMstr '_2_diff_FC2.pdf'],'ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(FIG2{3},['~/Downloads/' svstr '_' CMstr '_3_both_FC2.pdf'],'ContentType','vector');

end

% Visualization: deltaR2 mean summary for all FC
  for k=1:3
      opts.viz.usek=k;
      figs = bicsdual.viz_BICS_dual_allFC_summary(Delta_all, stats_all, ylbl, pinfo, opts.viz);
  end

  
% Visualization: elasticity
  figsE = bicsdual.viz_BICS_dual_elast(elasticity_all, stats_all, ylbl, pinfo, opts.vizE);
  set(figsE.elast_ribbons_Meq1,'Position',[-728 505 729 292]); title ''; xlabel ''; xtickangle(30);  % line
  set(figsE.elast_heatmaps_Meq1,'Position',[-320 54 321 743])                                        % heatmaps
  

% Visualization: surf & S-A axis corr
  opts.vizN.matchClim=1; 
  opts.vizN.surf_level='L1_route';
  % opts.vizN.surf_level='L1_diff';
  % opts.vizN.surf_level='L2_both';
  figsN = bicsdual.sa_corr_nodal(stats_all, SA_map, ylbl, pinfo, opts.vizN);
  set(figsN.fig_surf(1),'Position',[-1.7778   -0.0156    0.9000    0.9000])
  set(figsN.fig_surf(2),'Position',[-1.1444   -0.0156    0.9000    0.9000])
  set(figsN.fig_bars,'Position',[-2559 -539 729 245]); AX=findall(figsN.fig_bars,'type','axes');
  title(AX,''); xlabel(AX,''); ylabel(AX,''); xtickangle(AX,30);  % line



% Spin testof S-A axis corr
  iFC = 1;
  lbls_levels = {'L1_route' 'L1_diff' 'L2_both'};
  D=squeeze(stats_all.(lbls_levels{1}).node_deltaR2(1,:,iFC));
  [p_spin, r_obs, r_perm, info] = spin_pval_saaxis(D, SA_map, pinfo, 1, ...
    'n_perm', 100, 'stat', 'pearson', 'tail', 'two', 'seed', 42, 'method', 'vertex');


% --- Extract all global deltaR2 and plot bars
  ppostmp_bars=[-2559 594 643 203];
  label_level = {'L1_route' 'L1_diff' 'L2_both'}; useD=stats_all;
  DPLT=nan(3,Nfc); PVALS=nan(3,Nfc);
% Extract data
  for l = 1:3
      DPLT(l,:)=useD.(label_level{l}).global_deltaR2;
      PVALS(l,:)=useD.(label_level{l}).global_p;
  end
  DPLT(PVALS>0.05)=0;                                                       % filter by pvals
% Plot
  fig = myfig([CMstr ' (' intstr '): \DeltaR^2'],ppostmp_bars); hold on
  B = bar(DPLT', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  % for k=1:numel(B), set(B(k),'FaceColor',cols(k,:)); end
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz,usefont); % ylim([-.99 .99]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s','\DeltaR^2'));
  % legend(strrep(label_level,'_','-'), 'Location','northeastoutside'); 



%% =======      Specific bar plots from saved data      ======= %%
  pinfo.clabels_short   = LBL_ntwk;
  testset               = {[1 5] [1 6] [2 5] [2 6]};                        % SPE-CMY  SPE-DE  NE-CMY  NE-DE
  pair_idx              = cell2mat(cellfun(@(c) c(:),testset,'UniformOutput',0))';
  label_interact        = {'int_none' 'int_caliber' 'int_ed' 'int_both' 'total'};
  label_level           = {'L1_route' 'L1_diff' 'L2_both'};
  ppostmp_bars          = [-2559 594 588 203];
  ppostmp_bars_long     = [-2559 594 816 203];
  % tmptag='results_';
  tmptag='results_total_tmp';
  spinopt={pinfo,1,'n_perm',100,'stat','spearman','tail','two','seed',42,'method','vertex'};



% --------- SWITCH FOR INTERACTIONS
  ii=5; intstr=label_interact{ii};

  if ii<5, stats_all.meta.effect_mode='incremental';
  else,    stats_all.meta.effect_mode='total';
  end

% --- Global deltaR2, route & diff separate
  label_targ='global_deltaR2'; label_targ2='global_p';
  Dr=nan(2,Nfc); Dd=nan(2,Nfc); Pr=nan(2,Nfc); Pd=nan(2,Nfc);  
  it=0; iR=1; iD=2;
  for xx=[1 4]
          it = it + 1;
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-');                               % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          opts=S.meta.opts; label_comm_models=S.meta.label_comm_models; myelin_preds=S.meta.myelin_predictors{:};
          Delta_all=S.(intstr).Delta_all; P_all=S.(intstr).P_all; stats_all=S.(intstr).stats_all;
          SI_all=S.(intstr).SI_all; elasticity_all=S.(intstr).elasticity_all; pert_all=S.(intstr).pert_all;
        % Extract data
          Dr(it,:)=stats_all.(label_level{iR}).(label_targ);  Pr(it,:)=stats_all.(label_level{iR}).(label_targ2); % route      
          Dd(it,:)=stats_all.(label_level{iD}).(label_targ);  Pd(it,:)=stats_all.(label_level{iD}).(label_targ2); % diff
  end
% Check p-vals
  if any([Pr(:); Pd(:)]>=0.05), disp('Check pvals!!'); end 
% - Plot
% Routing
  fig = myfig(['Routing Models (' intstr '): \DeltaR^2'],ppostmp_bars); hold on
  B = bar(Dr', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  set(B(1),'FaceColor',cmap_6cm(1,:)); set(B(2),'FaceColor',cmap_6cm(2,:));  % SPE & NE
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); % ylim([-.99 .99]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s','\DeltaR^2')); 
  legend({'SPE' 'NE'}, 'Location','northeastoutside'); 
  set(gca,'TickLength',[0 0]); GCA1=gca;
% Diffusion
  fig = myfig(['Diffusion Models (' intstr '): \DeltaR^2'],ppostmp_bars); hold on
  B = bar(Dd', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  set(B(1),'FaceColor',cmap_6cm(5,:)); set(B(2),'FaceColor',cmap_6cm(6,:));  % CMY & DE
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); % ylim([-.99 .99]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s','\DeltaR^2')); 
  legend({'CMY' 'DE'}, 'Location','northeastoutside'); 
  set(gca,'TickLength',[0 0]); GCA2=gca;
  GCA1.Position=GCA2.Position;


  clear it iR iD Pr Pd Dr Dd


% --- All level 2 models (Route-Diff combined)
  label_targ='global_deltaR2'; label_targ2='global_p';
  Db=nan(4,Nfc); Pb=nan(4,Nfc); cmaptmp=gray(4+1);
  iB=3; PLTLBL=cell(1,4);
  for xx=1:4
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-'); PLTLBL{xx}=CMstr; % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          opts=S.meta.opts; label_comm_models=S.meta.label_comm_models; myelin_preds=S.meta.myelin_predictors{:};
          Delta_all=S.(intstr).Delta_all; P_all=S.(intstr).P_all; stats_all=S.(intstr).stats_all;
          SI_all=S.(intstr).SI_all; elasticity_all=S.(intstr).elasticity_all; pert_all=S.(intstr).pert_all;
        % Extract data
          Db(xx,:)=stats_all.(label_level{iB}).(label_targ);  Pb(xx,:)=stats_all.(label_level{iB}).(label_targ2); % route      
  end
% Check p-vals
  if any([Pb(:)]>=0.05), disp('Check pvals!!'); end
% Plot
  fig = myfig(['Rout-Diff Models (' intstr '): \DeltaR^2'],ppostmp_bars_long); hold on
  B = bar(Db', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  % for k=1:4, set(B(k),'FaceColor','k'); end
  for k=1:4, set(B(k),'FaceColor',cmaptmp(k,:)); end  
  % for k=1:4, B(k).LineWidth=.7; end
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); %ylim([0 0.2254]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s','\DeltaR^2')); 
  % [Hlow,Hup] = bar_dual_markers(B,Db,cmap_6cm,pair_idx,'barColors','k','marker','o','markerSize',100,'align','top');
  set(gca,'TickLength',[0 0]);

  clear iB Pb Db

% TMP SAVE
  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/1_int-none_1_route.pdf','ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/1_int-none_2_diff.pdf','ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/1_int-none_3_both.pdf','ContentType','vector');

  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/2_int-cal_1_route.pdf','ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/2_int-cal_2_diff.pdf','ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/2_int-cal_3_both.pdf','ContentType','vector');

  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/3_int-ed_1_route.pdf','ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/3_int-ed_2_diff.pdf','ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/3_int-ed_3_both.pdf','ContentType','vector');

  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/5_total_1_route.pdf','ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/5_total_2_diff.pdf','ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(gcf,'~/Downloads/5_total_3_both.pdf','ContentType','vector');


  % ylim([0 .26]); for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
  % ylim([0 0.079])
  % ylim([0 0.1300])



% --- Node-level: L1
  label_targ='node_deltaR2'; label_targ2='node_p';
  Dr=nan(2,Nnode,Nfc); Dd=nan(2,Nnode,Nfc); Pr=nan(2,Nnode,Nfc); Pd=nan(2,Nnode,Nfc);  
  it=0; iR=1; iD=2;
  for xx=[1 4]
          it = it + 1;
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-');                               % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          opts=S.meta.opts; label_comm_models=S.meta.label_comm_models; myelin_preds=S.meta.myelin_predictors{:};
          Delta_all=S.(intstr).Delta_all; P_all=S.(intstr).P_all; stats_all=S.(intstr).stats_all;
          SI_all=S.(intstr).SI_all; elasticity_all=S.(intstr).elasticity_all; pert_all=S.(intstr).pert_all;
        % Extract data
          Dr(it,:,:)=stats_all.(label_level{iR}).(label_targ);  Pr(it,:,:)=stats_all.(label_level{iR}).(label_targ2); % route      
          Dd(it,:,:)=stats_all.(label_level{iD}).(label_targ);  Pd(it,:,:)=stats_all.(label_level{iD}).(label_targ2); % diff
  end
% Check p-vals
  if any([Pr(:); Pd(:)]>=0.05), disp('Check pvals!!'); end 
  Dr(Pr>=0.05)=nan;   Dd(Pd>=0.05)=nan; 
% --- Corr with S-A axis
  Dr_corr=nan(2,Nfc); Dd_corr=nan(2,Nfc); Pr_corr=nan(2,Nfc); Pd_corr=nan(2,Nfc);
  for it = 1 :2 
      for ifc = 1 : Nfc
          Dt=squeeze(Dr(it,:,ifc)); [Pr_corr(it,ifc), Dr_corr(it,ifc),~,~]=spin_pval_saaxis(Dt,SA_map,spinopt{:}); % Route
          Dt=squeeze(Dd(it,:,ifc)); [Pd_corr(it,ifc), Dd_corr(it,ifc),~,~]=spin_pval_saaxis(Dt,SA_map,spinopt{:}); % Diff
      end 
  end
% - Plot
% Routing
  fig = myfig(['Routing Models (' intstr '): rho(SAaxis, deltaR2)'],ppostmp_bars); hold on
  B = bar(Dr_corr', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  set(B(1),'FaceColor',cmap_6cm(1,:)); set(B(2),'FaceColor',cmap_6cm(2,:));  % SPE & NE
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); ylim([-.75 .75]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s','\rho')); legend({'SPE' 'NE'}, 'Location','northeastoutside'); 
  set(gca,'TickLength',[0 0]);
% Diffusion
  fig = myfig(['Diffusion Models (' intstr '): rho(SAaxis, deltaR2)'],ppostmp_bars); hold on
  B = bar(Dd_corr', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  set(B(1),'FaceColor',cmap_6cm(5,:)); set(B(2),'FaceColor',cmap_6cm(6,:));  % CMY & DE
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); ylim([-.75 .75]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s','\rho')); legend({'CMY' 'DE'}, 'Location','northeastoutside'); 
  set(gca,'TickLength',[0 0]);
% -- Plot subset of surfaces
% % Routing
%   ifc=[1]; D=squeeze(Dr(1,:,ifc)); CMAX=ceil(prctile(D(:),95)*100)/100;
%   Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['SPE, deltaR2, ' intstr]);
%   setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);
% % Diffusion: CMY
%   ifc=[1]; D=squeeze(Dd(1,:,ifc)); CMAX=ceil(prctile(D(:),95)*100)/100;
%   Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['CMY, deltaR2, '  intstr]);
%   setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);
% % Diffusion: DE
%   ifc=[7 8]; D=squeeze(Dd(2,:,ifc)); CMAX=ceil(prctile(D(:),95)*100)/100;
%   Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['DE, deltaR2, '  intstr]);
%   setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);



% --- Node-level: L2
  label_targ='node_deltaR2'; label_targ2='node_p';
  Db=nan(4,Nnode,Nfc);  Pb=nan(4,Nnode,Nfc);  cmaptmp=gray(4+1);
  iB=3; PLTLBL=cell(1,4);
  for xx=1:4
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-'); PLTLBL{xx}=CMstr; % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          opts=S.meta.opts; label_comm_models=S.meta.label_comm_models; myelin_preds=S.meta.myelin_predictors{:};
          Delta_all=S.(intstr).Delta_all; P_all=S.(intstr).P_all; stats_all=S.(intstr).stats_all;
          SI_all=S.(intstr).SI_all; elasticity_all=S.(intstr).elasticity_all; pert_all=S.(intstr).pert_all;
        % Extract data
          Db(xx,:,:)=stats_all.(label_level{iB}).(label_targ);  
          Pb(xx,:,:)=stats_all.(label_level{iB}).(label_targ2);   % both      
  end
% Check p-vals
  if any(Pb(:)>=0.05), disp('Check pvals!!'); end 
  Db(Pb>=0.05)=nan; 
% --- Corr with S-A axis
  Db_corr=nan(4,Nfc); Pb_corr=nan(4,Nfc);
  for xx = 1 : 2 
      for ifc = 1 : Nfc
          Dt=squeeze(Db(xx,:,ifc)); [Pb_corr(xx,ifc), Db_corr(xx,ifc),~,~]=spin_pval_saaxis(Dt,SA_map,spinopt{:}); % Both
      end 
  end
% - Plot
% SA corr
  fig = myfig(['Rout-Diff Models (' intstr '): rho(SAaxis, deltaR2)'],ppostmp_bars_long); hold on
  B = bar(Db_corr', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  % for k=1:4, set(B(k),'FaceColor','k'); end
  for k=1:4, set(B(k),'FaceColor',cmaptmp(k,:)); end
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); ylim([-.75 .75]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s','\rho')); set(gca,'TickLength',[0 0]);
  % [Hlow,Hup] = bar_dual_markers(B,Db_corr,cmap_6cm,pair_idx,'barColors','k','marker','o','markerSize',50,'align','top');
  % legend([Hlow([1 3])  Hup(1:2)], {'SPE' 'NE' 'CMY' 'DE'});
% -- Plot subset of surfaces
%   icm=1; ifc=[1]; D=squeeze(Db(icm,:,ifc)); CMAX=ceil(prctile(D(:),95)*100)/100;
%   Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), [PLTLBL{icm} ', deltaR2, '   intstr]);
%   setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);


% Save spin test results
S=[];
S.route.Dr_corr=Dr_corr;  S.route.Pr_corr=Pr_corr;
S.diff.Dd_corr= Dd_corr;  S.diff.Pd_corr= Pd_corr;
S.both.Db_corr= Db_corr;  S.both.Pb_corr= Pb_corr;
S.opt=spinopt;
% Save results
  svdir=[locDeriv '/x_dual_BICS/grp_' parc];   
  if ~isfolder(svdir); mkdir(svdir); end
  save([svdir '/node_spintest_' intstr '.mat'],'S')



% % % % PVALS for int-both
% % % Pr_corr = [ 0.7183    0.5554    0.6613    0.0699    0.0010    0.1658    0.3327    0.9660; ...
% % %             0.9560    0.5325    0.2527    0.0549    0.0330    0.1978    0.1009    0.2468];
% % % 
% % % Pd_corr = [0.0430    0.0929    0.8561    0.0679    0.0150    0.0060    0.1429    0.1848; ...
% % %            0.4685    0.1429    0.7942    0.1728    0.1089    0.1459    0.0060    0.0010];
% % % 
% % % Pb_corr = [0.3516    0.6324    0.9780    0.1209    0.0120    0.0150    0.2817    0.5415; ...
% % %            0.7203    0.2178    0.8921    0.1648    0.1998    0.2068    0.7383    0.2967; ...
% % %            0.7832    0.3187    0.7822    0.1179    0.0040    0.0140    0.1738    0.2857; ...
% % %            0.5235    0.6244    0.6014    0.1059    0.1099    0.0799    0.7812    0.9331];



% ------- VERSION INCLUDING ONLY   SPE, CMY & DE
ppostmp=[-2559 594 701 203];
% --- Node-level: L1
  label_targ='node_deltaR2'; label_targ2='node_p';
  Dr=nan(2,Nnode,Nfc); Dd=nan(2,Nnode,Nfc); Pr=nan(2,Nnode,Nfc); Pd=nan(2,Nnode,Nfc);  
  it=0; iR=1; iD=2;
  for xx=[1 2]
          it = it + 1;
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-');                               % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          opts=S.meta.opts; label_comm_models=S.meta.label_comm_models; myelin_preds=S.meta.myelin_predictors{:};
          Delta_all=S.(intstr).Delta_all; P_all=S.(intstr).P_all; stats_all=S.(intstr).stats_all;
          SI_all=S.(intstr).SI_all; elasticity_all=S.(intstr).elasticity_all; pert_all=S.(intstr).pert_all;
        % Extract data
          Dr(it,:,:)=stats_all.(label_level{iR}).(label_targ);  Pr(it,:,:)=stats_all.(label_level{iR}).(label_targ2); % route      
          Dd(it,:,:)=stats_all.(label_level{iD}).(label_targ);  Pd(it,:,:)=stats_all.(label_level{iD}).(label_targ2); % diff
  end
% Check p-vals
  if any([Pr(:); Pd(:)]>=0.05), disp('Check pvals!!'); end 
  Dr(Pr>=0.05)=nan;   Dd(Pd>=0.05)=nan; 
  Dr(2,:,:)=[]; Pr(2,:,:)=[];
% --- Corr with S-A axis
  Dr_corr=nan(1,Nfc); Dd_corr=nan(2,Nfc); Pr_corr=nan(1,Nfc); Pd_corr=nan(2,Nfc);
% Diff
  for it = 1 :2 
      for ifc = 1 : Nfc
          Dt=squeeze(Dd(it,:,ifc)); [Pd_corr(it,ifc), Dd_corr(it,ifc),~,~]=spin_pval_saaxis(Dt,SA_map,spinopt{:}); % Diff
      end 
  end
% Route
  for ifc = 1 : Nfc
        Dt=squeeze(Dr(1,:,ifc)); [Pr_corr(1,ifc), Dr_corr(1,ifc),~,~]=spin_pval_saaxis(Dt,SA_map,spinopt{:}); % Route
  end
% - Plot
% Individual models (Routing & Diff)
  DPLT = [Dr_corr; Dd_corr];
  fig = myfig(['Individual Models (' intstr '): rho(SAaxis, deltaR2)'],ppostmp); hold on
  B = bar(DPLT', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  set(B(1),'FaceColor',cmap_6cm(1,:)); set(B(2),'FaceColor',cmap_6cm(5,:)); set(B(3),'FaceColor',cmap_6cm(6,:));  % SPE, CMY, DE
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); ylim([-.75 .75]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s','\rho')); legend({'SPE' 'CMY' 'DE'}, 'Location','northeastoutside'); 
  set(gca,'TickLength',[0 0]);
% -- Plot subset of surfaces: hetero clims
  ifc=[8]; 
% Routing
  % ifc=[1]; 
  D=squeeze(Dr(1,:,ifc)); CMAX=ceil(prctile(D(:),95)*100)/100;
  Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['SPE, deltaR2, ' intstr]);
  setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);
% Diffusion: CMY
  % ifc=[1]; 
  D=squeeze(Dd(1,:,ifc)); CMAX=ceil(prctile(D(:),95)*100)/100;
  Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['CMY, deltaR2, '  intstr]);
  setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);
% Diffusion: DE
  % ifc=[1]; 
  D=squeeze(Dd(2,:,ifc)); CMAX=ceil(prctile(D(:),95)*100)/100;
  Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['DE, deltaR2, '  intstr]);
  setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);

% -- Plot subset of surfaces: common clim
  % ifc=[1]; CMAX=0.35;
  % ifc=[5]; CMAX=0.35;
  ifc=[8]; CMAX=0.35;
  % ifc=[1 5 8]; CMAX=0.35;
  svstr=[num2str(ifc) '-' strrep(ylbl{ifc},'_','') '_'];
% Routing
  D=squeeze(Dr(1,:,ifc)); 
  Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['SPE, deltaR2, ' intstr]);
  setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr 'SPE.pdf'],'ContentType','vector');
% Diffusion: CMY
  D=squeeze(Dd(1,:,ifc)); 
  Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['CMY, deltaR2, '  intstr]);
  setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr 'CMY.pdf'],'ContentType','vector');
% Diffusion: DE
  D=squeeze(Dd(2,:,ifc));
  Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['DE, deltaR2, '  intstr]);
  setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr 'DE.pdf'],'ContentType','vector');

  

% --- Node-level: L2
  label_targ='node_deltaR2'; label_targ2='node_p';
  Db=nan(2,Nnode,Nfc);  Pb=nan(2,Nnode,Nfc);  cmaptmp=gray(2+1);
  iB=3; PLTLBL=cell(1,4);
  for xx=1:2
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-'); PLTLBL{xx}=CMstr; % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          opts=S.meta.opts; label_comm_models=S.meta.label_comm_models; myelin_preds=S.meta.myelin_predictors{:};
          Delta_all=S.(intstr).Delta_all; P_all=S.(intstr).P_all; stats_all=S.(intstr).stats_all;
          SI_all=S.(intstr).SI_all; elasticity_all=S.(intstr).elasticity_all; pert_all=S.(intstr).pert_all;
        % Extract data
          Db(xx,:,:)=stats_all.(label_level{iB}).(label_targ);  
          Pb(xx,:,:)=stats_all.(label_level{iB}).(label_targ2);   % both      
  end
% Check p-vals
  if any(Pb(:)>=0.05), disp('Check pvals!!'); end 
  Db(Pb>=0.05)=nan; 
% --- Corr with S-A axis
  Db_corr=nan(2,Nfc); Pb_corr=nan(2,Nfc);
  for xx = 1 : 2 
      for ifc = 1 : Nfc
          Dt=squeeze(Db(xx,:,ifc)); [Pb_corr(xx,ifc), Db_corr(xx,ifc),~,~]=spin_pval_saaxis(Dt,SA_map,spinopt{:}); % Both
      end 
  end
% - Plot
% SA corr
  fig = myfig(['Rout-Diff Models (' intstr '): rho(SAaxis, deltaR2)'],ppostmp); hold on
  B = bar(Db_corr', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  % for k=1:4, set(B(k),'FaceColor','k'); end
  for k=1:2, set(B(k),'FaceColor',cmaptmp(k,:)); end
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); ylim([-.75 .75]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s','\rho')); set(gca,'TickLength',[0 0]);
  legend({'SPE-CMY' 'SPE-DE'}, 'Location','northeastoutside');

% -- Plot subset of surfaces: hetero clims
  ifc=[8]; 
% SPE-CMY
  D=squeeze(Db(1,:,ifc)); CMAX=ceil(prctile(D(:),95)*100)/100;
  Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['SPE-CMY, deltaR2, ' intstr]);
  setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);
% SPE-DE
  D=squeeze(Db(2,:,ifc)); CMAX=ceil(prctile(D(:),95)*100)/100;
  Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['SPE-DE, deltaR2, ' intstr]);
  setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);

% -- Plot subset of surfaces: common clims
  % ifc=[1]; CMAX=0.35;
  % ifc=[5]; CMAX=0.35;
  ifc=[8]; CMAX=0.35;
  % ifc=[1 5 8]; CMAX=0.35;
  svstr=[num2str(ifc) '-' strrep(ylbl{ifc},'_','') '_'];
% SPE-CMY
  D=squeeze(Db(1,:,ifc)); 
  Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['SPE-CMY, deltaR2, ' intstr]);
  setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr 'SPE-CMY.pdf'],'ContentType','vector');
% SPE-DE
  D=squeeze(Db(2,:,ifc)); 
  Hcx = plot_conn_surf(D', pinfo, 'cortex', ylbl(ifc), ['SPE-DE, deltaR2, ' intstr]);
  setsurf(Hcx, [0 CMAX], [.7 .7 .7; cm4]);
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr 'SPE-DE.pdf'],'ContentType','vector');




% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % %% =======   LOOP VERSION OF Modeling 1c: DUAL Blockwise RSN FC prediction
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % All tests designed to include 1 routing & 1 diffusion model together
% % % % Interactions added as well
% % % 
% % % % --- Setup
% % % testset={[1 5] [1 6] [2 5] [ 2 6]};
% % % for xx=1:4
% % %       CC=testset{xx};
% % %       CMstr                 = strjoin(CMLBL_simple(CC),'-');
% % %       ALLY                  = Yz_m;
% % %       EDz                   = nzzscore(EDnrm,0);
% % %       ALLX                  = [ALLCM {ALLCM_delay} {ALLCM_bin} {EDz}]; 
% % %       pinfo.clabels_short   = LBL_ntwk;
% % %       uselower              = 1; 
% % %       elast_mode            = 'zrow';
% % %       S=[];
% % % 
% % % 
% % %     % Run Opts
% % %       opts = struct( ...
% % %                 'stageA', struct( ...
% % %                             'use_lower_only', uselower, ...
% % %                             'alpha', 0.05 ...                                   % intended for FDR (add later)
% % %                                 ), ...
% % %                 'stageB', struct( ...
% % %                             'nperm', 200, ...
% % %                             'eps_scale', 0.05, ...
% % %                             'use_lower_only', uselower, ...
% % %                             'rng_seed', 13, ...
% % %                             'perm_mode', 'yresid', ...
% % %                             'elasticity_mode', elast_mode ...                   % 'xscale' | 'zadd'
% % %                             ), ...
% % %                 'viz', struct( ...                                              % main results
% % %                             'figpos1', [-2559 54 563 743], ...                  % heatmaps (M=1 case)
% % %                             'fc_agg', 'mean', ...                               % 'mean'|'median'
% % %                             'labelON', true, ...                                % turn some plot labels on/off
% % %                             'elasticity', [], ...                               % for Figure 3 style panels
% % %                             'elasticity_norm', 'none', ...                      % 'none'|'by_delta'
% % %                             'elasticity_norm_floor', 1e-6, ...
% % %                             'elasticity_as_percent', true, ...
% % %                             'ribbons_xdim', 'bands' ...                         % 'models' | 'bands'
% % %                             ), ...
% % %                 'vizE', struct( ...                                             % elastcitiy results
% % %                             'figpos1', [-2559 54 321 743], ...                  % heatmaps (M=1 case)
% % %                             'fc_agg', 'mean', ...                               % 'mean'|'median'
% % %                             'labelON', true, ...                                % turn some plot labels on/off
% % %                             'elasticity_mode', elast_mode, ...                  % 'xscale' | 'zadd'
% % %                             'norm', 'none', ...                                 % 'none'|'by_delta'
% % %                             'norm_floor', 1e-6, ...
% % %                             'ribbons_xdim', 'bands' ...                         % 'models' | 'bands'
% % %                             ), ...
% % %                 'vizN', struct( ...                                             % node results
% % %                             'figpos_bars', [-936 552 936 245], ...              % M=1 case
% % %                             'corr_type', 'Spearman', ...                        % 'mean'|'median'
% % %                             'absval', false, ...                                % 'xscale' | 'zadd'
% % %                             'm_eq1_colors', [], ...                             % 'none'|'by_delta'
% % %                             'bands', 1:Nfc, ...                                 % numeric
% % %                             'labelON', true, ...                                % turn some plot labels on/off
% % %                             'make_surfaces', true, ...
% % %                             'surf_bands', [], ...
% % %                             'surf_level', 'L2_both',...
% % %                             'surf_title', 'deltaR2', ...
% % %                             'surf_cmax', [], ...                                % auto (95th pct of |ΔR²|)
% % %                             'surf_cmap', [.7 .7 .7; cm4] ... 
% % %                             ) ...
% % %                     );
% % % 
% % %     % -------------------------------------------------------------------------
% % %     % =====  Base (no interactions) =====
% % %     % -------------------------------------------------------------------------
% % % 
% % %     tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
% % %                    'label_comm_models',CMstr, ...
% % %                    'modes', 'all', 'levels', [1 2],                   ...
% % %                    'interaction', 'none');
% % % 
% % % 
% % %     % STAGE A: Compute ΔR² per RSN block, model spec, and FC target.
% % %       [Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);
% % % 
% % %     % STAGE B: Sensitivity and elasticity testing
% % %       [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts.stageB);
% % % 
% % %     % Save it
% % %       S.int_none.Delta_all=Delta_all;
% % %       S.int_none.P_all=P_all;
% % %       S.int_none.stats_all=stats_all;
% % %       S.int_none.SI_all=SI_all;
% % %       S.int_none.elasticity_all=elasticity_all;
% % %       S.int_none.pert_all=pert_all;
% % % 
% % % 
% % %     % -------------------------------------------------------------------------
% % %     % ==========          CALIBER interactions  (L1 and L2)          ==========
% % %     % -------------------------------------------------------------------------
% % % 
% % %     tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
% % %                    'label_comm_models',CMstr, ...
% % %                    'modes', 'all', 'levels', [1 2],                   ...
% % %                    'interaction', 'caliber', 'interact_at', 'both');
% % % 
% % %     % STAGE A: Compute ΔR² per RSN block, model spec, and FC target.
% % %       [Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);
% % % 
% % %     % STAGE B: Sensitivity and elasticity testing
% % %       [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts.stageB);
% % % 
% % %     % Save it
% % %       S.int_caliber.Delta_all=Delta_all;
% % %       S.int_caliber.P_all=P_all;
% % %       S.int_caliber.stats_all=stats_all;
% % %       S.int_caliber.SI_all=SI_all;
% % %       S.int_caliber.elasticity_all=elasticity_all;
% % %       S.int_caliber.pert_all=pert_all;
% % % 
% % %     % -------------------------------------------------------------------------
% % %     % ==========             ED interactions  (L1 and L2)            ==========
% % %     % -------------------------------------------------------------------------
% % % 
% % %     tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
% % %                    'label_comm_models',CMstr, ...
% % %                    'modes', 'all', 'levels', [1 2],                   ...
% % %                    'interaction', 'ed', 'interact_at', 'both');
% % % 
% % %     % STAGE A: Compute ΔR² per RSN block, model spec, and FC target.
% % %       [Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);
% % % 
% % %     % STAGE B: Sensitivity and elasticity testing
% % %       [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts.stageB);
% % % 
% % %     % Save it
% % %       S.int_ed.Delta_all=Delta_all;
% % %       S.int_ed.P_all=P_all;
% % %       S.int_ed.stats_all=stats_all;
% % %       S.int_ed.SI_all=SI_all;
% % %       S.int_ed.elasticity_all=elasticity_all;
% % %       S.int_ed.pert_all=pert_all;
% % % 
% % %     % -------------------------------------------------------------------------
% % %     % ==========            BOTH interactions  (L1 and L2)           ==========
% % %     % -------------------------------------------------------------------------
% % % 
% % %     tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
% % %                    'label_comm_models',CMstr, ...
% % %                    'modes', 'all', 'levels', [1 2],                   ...
% % %                    'interaction', 'both', 'interact_at', 'both');
% % % 
% % %     % STAGE A: Compute ΔR² per RSN block, model spec, and FC target.
% % %       [Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);
% % % 
% % %     % STAGE B: Sensitivity and elasticity testing
% % %       [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts.stageB);
% % % 
% % % 
% % %     % Save it
% % %       S.int_both.Delta_all=Delta_all;
% % %       S.int_both.P_all=P_all;
% % %       S.int_both.stats_all=stats_all;
% % %       S.int_both.SI_all=SI_all;
% % %       S.int_both.elasticity_all=elasticity_all;
% % %       S.int_both.pert_all=pert_all;
% % % 
% % %       S.meta.opts=opts;
% % %       S.meta.label_comm_models=strjoin(CMLBL_simple(CC),'-');
% % %       S.meta.myelin_predictors={{'MTsat','gratio','delay'}};
% % % 
% % %     % Save results
% % %       svdir=[locDeriv '/x_dual_BICS/grp_nos_' parc];   
% % %       if ~isfolder(svdir); mkdir(svdir); end
% % %       save([svdir '/results_' CMstr '.mat'],'S')
% % % 
% % % end
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % 
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % %% =======   LOOP VERSION 2 OF Modeling 1c: DUAL Blockwise RSN FC prediction
% % % %% =======   ADDS TOTAL CONTRIBUTION CASE
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % All tests designed to include 1 routing & 1 diffusion model together
% % % % Interactions added as well
% % % 
% % % % --- Setup
% % % testset={[1 5] [1 6] [2 5] [ 2 6]};
% % % svdir=[locDeriv '/x_dual_BICS/grp_nos_' parc];   
% % % if ~isfolder(svdir); mkdir(svdir); end
% % % 
% % % for xx=1:4
% % %       CC=testset{xx};
% % %       CMstr                 = strjoin(CMLBL_simple(CC),'-');
% % %       ALLY                  = Yz_m;
% % %       EDz                   = nzzscore(EDnrm,0);
% % %       ALLX                  = [ALLCM {ALLCM_delay} {ALLCM_bin} {EDz}]; 
% % %       pinfo.clabels_short   = LBL_ntwk;
% % %       uselower              = 1; 
% % %       elast_mode            = 'zrow';
% % % 
% % % 
% % %       load([svdir '/results_' CMstr '.mat'])
% % % 
% % % 
% % %     % Run Opts
% % %       opts = struct( ...
% % %                 'stageA', struct( ...
% % %                             'use_lower_only', uselower, ...
% % %                             'alpha', 0.05 ...                                   % intended for FDR (add later)
% % %                                 ), ...
% % %                 'stageB', struct( ...
% % %                             'nperm', 200, ...
% % %                             'eps_scale', 0.05, ...
% % %                             'use_lower_only', uselower, ...
% % %                             'rng_seed', 13, ...
% % %                             'perm_mode', 'yresid', ...
% % %                             'elasticity_mode', elast_mode ...                   % 'xscale' | 'zadd'
% % %                             ), ...
% % %                 'viz', struct( ...                                              % main results
% % %                             'figpos1', [-2559 54 563 743], ...                  % heatmaps (M=1 case)
% % %                             'fc_agg', 'mean', ...                               % 'mean'|'median'
% % %                             'labelON', true, ...                                % turn some plot labels on/off
% % %                             'elasticity', [], ...                               % for Figure 3 style panels
% % %                             'elasticity_norm', 'none', ...                      % 'none'|'by_delta'
% % %                             'elasticity_norm_floor', 1e-6, ...
% % %                             'elasticity_as_percent', true, ...
% % %                             'ribbons_xdim', 'bands' ...                         % 'models' | 'bands'
% % %                             ), ...
% % %                 'vizE', struct( ...                                             % elastcitiy results
% % %                             'figpos1', [-2559 54 321 743], ...                  % heatmaps (M=1 case)
% % %                             'fc_agg', 'mean', ...                               % 'mean'|'median'
% % %                             'labelON', true, ...                                % turn some plot labels on/off
% % %                             'elasticity_mode', elast_mode, ...                  % 'xscale' | 'zadd'
% % %                             'norm', 'none', ...                                 % 'none'|'by_delta'
% % %                             'norm_floor', 1e-6, ...
% % %                             'ribbons_xdim', 'bands' ...                         % 'models' | 'bands'
% % %                             ), ...
% % %                 'vizN', struct( ...                                             % node results
% % %                             'figpos_bars', [-936 552 936 245], ...              % M=1 case
% % %                             'corr_type', 'Spearman', ...                        % 'mean'|'median'
% % %                             'absval', false, ...                                % 'xscale' | 'zadd'
% % %                             'm_eq1_colors', [], ...                             % 'none'|'by_delta'
% % %                             'bands', 1:Nfc, ...                                 % numeric
% % %                             'labelON', true, ...                                % turn some plot labels on/off
% % %                             'make_surfaces', true, ...
% % %                             'surf_bands', [], ...
% % %                             'surf_level', 'L2_both',...
% % %                             'surf_title', 'deltaR2', ...
% % %                             'surf_cmax', [], ...                                % auto (95th pct of |ΔR²|)
% % %                             'surf_cmap', [.7 .7 .7; cm4] ... 
% % %                             ) ...
% % %                     );
% % % 
% % % 
% % %     % -------------------------------------------------------------------------
% % %     % ==========            TOTAL CONTRIBUTION           ==========
% % %     % -------------------------------------------------------------------------
% % % 
% % %     tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
% % %                    'label_comm_models',CMstr, ...
% % %                    'modes', 'all', 'levels', [1 2],                   ...
% % %                    'interaction', 'both', 'interact_at', 'both', 'effect_mode', 'total');
% % % 
% % %     % STAGE A: Compute ΔR² per RSN block, model spec, and FC target.
% % %       tic;
% % %       [Delta_all, P_all, stats_all] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);
% % %       disp(['End stageA: ' CMstr ' in ' num2str(toc) 's'])
% % % 
% % %     % STAGE B: Sensitivity and elasticity testing
% % %       tic;
% % %       [SI_all, elasticity_all, pert_all] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts.stageB);
% % %       disp(['End stageB: ' CMstr ' in ' num2str(toc) 's'])
% % % 
% % % 
% % % 
% % %     % Save it
% % %       S.total.Delta_all=Delta_all;
% % %       S.total.P_all=P_all;
% % %       S.total.stats_all=stats_all;
% % %       S.total.SI_all=SI_all;
% % %       S.total.elasticity_all=elasticity_all;
% % %       S.total.pert_all=pert_all;
% % % 
% % %       S.meta.opts=opts;
% % %       S.meta.label_comm_models=strjoin(CMLBL_simple(CC),'-');
% % %       S.meta.myelin_predictors={{'MTsat','gratio','delay'}};
% % % 
% % %     % Save results
% % %       save([svdir '/results_total_tmp' CMstr '.mat'],'S')
% % % 
% % % end








% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%% =======   Modeling 1d: Myelin vs Caliber Dual BICS
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% All tests designed to include 1 routing & 1 diffusion model together

% --- Setup
  CC                    = [1 5];
  CMstr                 = strjoin(CMLBL_simple(CC),'-');
  ALLY                  = Yz_m;
  EDz                   = nzzscore(EDnrm,0);
  ALLX                  = [ALLCM {ALLCM_delay} {ALLCM_bin} {EDz}]; 
  pinfo.clabels_short   = LBL_ntwk;
  uselower              = 1; 
  elast_mode            = 'zrow';
  FCAGG                 = 'median';

% DEBUGGING
  % % % optsA.debug = struct('probe',true,'probe_blocks',false,'probe_nodes',false,'trace_menus',true,'save','dbg_stageA.mat');
  % % % optsB.debug = struct('probe_B',true,'save','dbg_stageB.mat');
  % % % tests = struct('myelin_predictors', {{'MTsat','gratio','delay'}}, ...
  % % %              'label_comm_models',strjoin(CMLBL_simple(CC),'-'), ...
  % % %              'modes', 'all', 'levels', [1 2],                   ...
  % % %              'interaction', 'caliber', 'interact_at', 'both');
  % % % 
  % % %              % 'interaction', 'none');
  % % % [Delta_all,~,~] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, optsA);
  % % % [~,~,~] = bicsdual.fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, optsB);

% Run Opts
  opts = struct( ...
            'stageA', struct( ...
                        'use_lower_only', uselower, ...
                        'alpha', 0.05 ...                                   % intended for FDR (add later)
                            ), ...
            'stageB', struct( ...
                        'nperm', 200, ...
                        'eps_scale', 0.05, ...
                        'use_lower_only', uselower, ...
                        'rng_seed', 13, ...
                        'perm_mode', 'yresid', ...
                        'elasticity_mode', elast_mode ...                   % 'xscale' | 'zadd'
                        ), ...
            'viz', struct( ...                                              % main results
                        'figpos1', [-2559 54 563 743], ...                  % heatmaps (M=1 case)
                        'labelON', true, ...                                % turn some plot labels on/off
                        'fc_agg', FCAGG, ...                                % 'mean'|'median'
                        'elasticity', [], ...                               % for Figure 3 style panels
                        'elasticity_norm', 'none', ...                      % 'none'|'by_delta'
                        'elasticity_norm_floor', 1e-6, ...
                        'elasticity_as_percent', true, ...
                        'ribbons_xdim', 'bands' ...                         % 'models' | 'bands'
                        ), ...
            'vizE', struct( ...                                             % elastcitiy results
                        'figpos1', [-2559 54 321 743], ...                  % heatmaps (M=1 case)
                        'fc_agg', FCAGG, ...                                % 'mean'|'median'
                        'labelON', true, ...                                % turn some plot labels on/off
                        'elasticity_mode', elast_mode, ...                  % 'xscale' | 'zadd'
                        'norm', 'none', ...                                 % 'none'|'by_delta'
                        'norm_floor', 1e-6, ...
                        'ribbons_xdim', 'bands' ...                         % 'models' | 'bands'
                        ), ...
            'vizN', struct( ...                                             % node results
                        'figpos_bars', [-936 552 936 245], ...              % M=1 case
                        'corr_type', 'Spearman', ...                        % 'mean'|'median'
                        'absval', false, ...                                % 'xscale' | 'zadd'
                        'm_eq1_colors', [], ...                             % 'none'|'by_delta'
                        'bands', 1:Nfc, ...                                 % numeric
                        'labelON', true, ...                                % turn some plot labels on/off
                        'make_surfaces', true, ...
                        'surf_bands', [], ...
                        'surf_level', 'L2_both',...
                        'surf_title', 'deltaR2', ...
                        'surf_cmax', [], ...                                % auto (95th pct of |ΔR²|)
                        'surf_cmap', [.7 .7 .7; cm4] ... 
                        ) ...
                );

% -------------------------------------------------------------------------
% =====  Test 1:  Myelin vs base (binary_cm + ED)
% -------------------------------------------------------------------------

use_myelin_str = 'mtsat';

tests = struct('myelin_predictors', {{use_myelin_str}}, 'base_predictors', {{'binary'}}, ...
               'label_comm_models',CMstr, ...
               'modes', 'single', 'levels', [1 2],                   ...
               'interaction', 'none', 'effect_mode', 'incremental');


% Test 1: Myelin vs base
  [Delta_all_myelin,  P_all_myelin,  stats_all_myelin] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);

% Test 2: Caliber vs base
  tests.myelin_predictors = {'caliber'};
  [Delta_all_caliber, P_all_caliber, stats_all_caliber] = bicsdual.fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts.stageA);
    
% --- Compute delta deltaR2
  fields = {'L1_route' 'L1_diff' 'L2_both'}; 
  delta_stats_all = stats_all_myelin;                   % use myelin stats structure as start point
  for k = 1:3
      % RSN
        delta_deltaR2.(fields{k})   = Delta_all_myelin.(fields{k}) - Delta_all_caliber.(fields{k});
      % -- Global
      % deltaR2
        tmpfield = 'global_deltaR2';
        x = stats_all_myelin.(fields{k}).(tmpfield);
        y = stats_all_caliber.(fields{k}).(tmpfield);
        delta_stats_all.(fields{k}).(tmpfield) = x - y;                     % difference
      % pvals
        tmpfield = 'global_p';
        x = stats_all_myelin.(fields{k}).(tmpfield);
        y = stats_all_caliber.(fields{k}).(tmpfield);
        delta_stats_all.(fields{k}).(tmpfield) = max(x, y);               % max
      % -- ntwk
      % deltaR2
        tmpfield = 'ntwk_deltaR2';
        x = stats_all_myelin.(fields{k}).(tmpfield);
        y = stats_all_caliber.(fields{k}).(tmpfield);
        delta_stats_all.(fields{k}).(tmpfield) = x - y;                     % difference
      % pvals
        tmpfield = 'ntwk_p';
        x = stats_all_myelin.(fields{k}).(tmpfield);
        y = stats_all_caliber.(fields{k}).(tmpfield);
        delta_stats_all.(fields{k}).(tmpfield) = max(x, y);                 % max
      % -- node
      % deltaR2
        tmpfield = 'node_deltaR2';
        x = stats_all_myelin.(fields{k}).(tmpfield);
        y = stats_all_caliber.(fields{k}).(tmpfield);
        delta_stats_all.(fields{k}).(tmpfield) = x - y;                     % difference
      % pvals
        tmpfield = 'node_p';
        x = stats_all_myelin.(fields{k}).(tmpfield);
        y = stats_all_caliber.(fields{k}).(tmpfield);
        delta_stats_all.(fields{k}).(tmpfield) = max(x, y);                 % max
  end
  delta_stats_all.meta.test = [use_myelin_str ' vs caliber'];
  S=[];
  S.delta_deltaR2 = delta_deltaR2;
  S.delta_stats_all = delta_stats_all;

% Save results
  svdir=[locDeriv '/x_dual_BICS/grp_' parc];   
  save([svdir '/1_' use_myelin_str '-vs-caliber_' CMstr '.mat'],'S')


%% =======      Specific bar plots from saved data      ======= %%
  pinfo.clabels_short   = LBL_ntwk;
  testset               = {[1 5] [1 6] [2 5] [2 6]};                        % SPE-CMY  SPE-DE  NE-CMY  NE-DE
  pair_idx              = cell2mat(cellfun(@(c) c(:),testset,'UniformOutput',0))';
  label_interact        = {'int_none' 'int_caliber' 'int_ed' 'int_both' 'total'};
  label_level           = {'L1_route' 'L1_diff' 'L2_both'};
  ppostmp_bars          = [-2559 594 588 203];
  ppostmp_bars_long     = [-2559 594 816 203];

% SWITCH FOR MYELIN DATA
  use_myelin_str = 'mtsat';
  tmptag = ['1_' use_myelin_str '-vs-caliber_'];


% ------ Global

% --- Level 1: Route & Diff separate
  label_targ='global_deltaR2'; label_targ2='global_p'; tmpstr=[use_myelin_str ' - Caliber'];
  Dr=nan(2,Nfc); Dd=nan(2,Nfc); Pr=nan(2,Nfc); Pd=nan(2,Nfc);  
  it=0; iR=1; iD=2;
  for xx=[1 4]
          it = it + 1;
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-');                               % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          label_comm_models=S.delta_stats_all.meta.label_comm_models;  % opts=S.meta.opts;
          deltaR2   = S.delta_deltaR2; 
          stats_all = S.delta_stats_all;
        % Extract data
          Dr(it,:)=stats_all.(label_level{iR}).(label_targ);  Pr(it,:)=stats_all.(label_level{iR}).(label_targ2); % route      
          Dd(it,:)=stats_all.(label_level{iD}).(label_targ);  Pd(it,:)=stats_all.(label_level{iD}).(label_targ2); % diff
  end
% Check p-vals
  if any([Pr(:); Pd(:)]>=0.05), disp('Check pvals!!'); end 
% - Plot
% Routing
  fig1 = myfig(['Routing Models (' tmpstr '): \Delta \DeltaR^2'],ppostmp_bars); hold on
  B = bar(Dr', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  set(B(1),'FaceColor',cmap_6cm(1,:)); set(B(2),'FaceColor',cmap_6cm(2,:));  % SPE & NE
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); % ylim([-.99 .99]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s%s','\Delta', '\DeltaR^2')); 
  legend({'SPE' 'NE'}, 'Location','northeastoutside'); 
  set(gca,'TickLength',[0 0]); GCA1=gca;
  ylim([-.2 .15]); for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end; legend off

% Diffusion
  fig2 = myfig(['Diffusion Models (' tmpstr '): \Delta \DeltaR^2'],ppostmp_bars); hold on
  B = bar(Dd', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  set(B(1),'FaceColor',cmap_6cm(5,:)); set(B(2),'FaceColor',cmap_6cm(6,:));  % CMY & DE
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); % ylim([-.99 .99]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s%s','\Delta', '\DeltaR^2')); 
  legend({'CMY' 'DE'}, 'Location','northeastoutside'); 
  set(gca,'TickLength',[0 0]); GCA2=gca;
  GCA1.Position=GCA2.Position;
  ylim([-.2 .15]); for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end; legend off

% Save
  svstr='1_global_';
  set(gcf,'Renderer','painters'); exportgraphics(fig1,['~/Downloads/' svstr tmptag '1-route.pdf'],'ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(fig2,['~/Downloads/' svstr tmptag '2-diff.pdf'],'ContentType','vector');




% --- Level 2: (Route-Diff combined)
  label_targ='global_deltaR2'; label_targ2='global_p'; tmpstr=[use_myelin_str ' - Caliber'];
  Db=nan(4,Nfc); Pb=nan(4,Nfc); cmaptmp=gray(4+1);
  iB=3; PLTLBL=cell(1,4);
  for xx=1:4
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-'); PLTLBL{xx}=CMstr; % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          label_comm_models=S.delta_stats_all.meta.label_comm_models;  % opts=S.meta.opts;
          deltaR2   = S.delta_deltaR2; 
          stats_all = S.delta_stats_all;
        % Extract data
          Db(xx,:)=stats_all.(label_level{iB}).(label_targ);  Pb(xx,:)=stats_all.(label_level{iB}).(label_targ2); % route      
  end
% Check p-vals
  if any([Pb(:)]>=0.05), disp('Check pvals!!'); end
% Plot
  fig3 = myfig(['Rout-Diff Models (' tmpstr '): \Delta \DeltaR^2'],ppostmp_bars_long); hold on
  B = bar(Db', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  % for k=1:4, set(B(k),'FaceColor','k'); end
  for k=1:4, set(B(k),'FaceColor',cmaptmp(k,:)); end  
  % for k=1:4, B(k).LineWidth=.7; end
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fntsz-2,usefont); %ylim([0 0.2254]);
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('%s%s','\Delta', '\DeltaR^2'));  
  % [Hlow,Hup] = bar_dual_markers(B,Db,cmap_6cm,pair_idx,'barColors','k','marker','o','markerSize',100,'align','top');
  set(gca,'TickLength',[0 0]);

% Save
  svstr='1_global_';
  set(gcf,'Renderer','painters'); exportgraphics(fig3,['~/Downloads/' svstr tmptag '3-both.pdf'],'ContentType','vector');



% ------ RSN level

% --- Level 1: Route & Diff separate
  label_targ='ntwk_deltaR2'; label_targ2='ntwk_p'; tmpstr=[use_myelin_str ' - Caliber'];
  Dr=nan(2,Nntwk,Nntwk,Nfc); Dd=nan(2,Nntwk,Nntwk,Nfc); 
  Pr=nan(2,Nntwk,Nntwk,Nfc); Pd=nan(2,Nntwk,Nntwk,Nfc);  
  it=0; iR=1; iD=2;
  for xx=[1 4]
          it = it + 1;
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-');                               % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          label_comm_models=S.delta_stats_all.meta.label_comm_models;  % opts=S.meta.opts;
          deltaR2   = S.delta_deltaR2; 
          stats_all = S.delta_stats_all;
        % Extract data
          Dr(it,:,:,:)=stats_all.(label_level{iR}).(label_targ);  Pr(it,:,:,:)=stats_all.(label_level{iR}).(label_targ2); % route      
          Dd(it,:,:,:)=stats_all.(label_level{iD}).(label_targ);  Pd(it,:,:,:)=stats_all.(label_level{iD}).(label_targ2); % diff
  end
% Check p-vals
  if any([Pr(:); Pd(:)]>=0.05), disp('Check pvals!!'); end 
% - Plot
  CMAX=0.2; %CMAX = prctile(abs([Dr(:); Dd(:)]),97);
% Routing
  CM_labels = {'SPE'  'NE'};
  for xx=1:2
      fig = myfig([CM_labels{xx} ' (' tmpstr '): \Delta \DeltaR^2'],[-2559 54 563 300+260*4]);
      tl = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
      for k=1:Nfc
          Dt = squeeze(Dr(xx,:,:,k)); ax = nexttile(tl,k);
          imagesc(ax, Dt); axis(ax,'square'); box(ax,'on'); title(ylbl{k});
          add_mat_box(ax, Dt, 'Color','k','LineWidth',0.5,'Outer',true);
          clim(ax, [-CMAX CMAX]); colormap(ax,cm2); cb=colorbar(ax); font(fntsz-4,usefont);
          set(ax,'YTick',1:Nntwk,'YTickLabels',LBL_ntwk,'XTick',1:Nntwk,'XTickLabels',LBL_ntwk); xtickangle(40); ytickangle(0);
      end
  end
% Diffusion
  CM_labels = {'CMY'  'DE'};
  for xx=1:2
      fig = myfig([CM_labels{xx} ' (' tmpstr '): \Delta \DeltaR^2'],[-2559 54 563 300+260*4]);
      tl = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
      for k=1:Nfc
          Dt = squeeze(Dd(xx,:,:,k)); ax = nexttile(tl,k);
          imagesc(ax, Dt); axis(ax,'square'); box(ax,'on'); title(ylbl{k});
          add_mat_box(ax, Dt, 'Color','k','LineWidth',0.5,'Outer',true);
          clim(ax, [-CMAX CMAX]); colormap(ax,cm2); cb=colorbar(ax); font(fntsz-4,usefont);
          set(ax,'YTick',1:Nntwk,'YTickLabels',LBL_ntwk,'XTick',1:Nntwk,'XTickLabels',LBL_ntwk); xtickangle(40); ytickangle(0);
      end
  end

  svstr='2_ntwk_';
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr tmptag '1-spe.pdf'],'ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr tmptag '3-cmy.pdf'],'ContentType','vector');
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr tmptag '4-de.pdf'],'ContentType','vector');




% --- Level 2: (Route-Diff combined)
  label_targ='ntwk_deltaR2'; label_targ2='ntwk_p'; tmpstr=[use_myelin_str ' - Caliber'];
  Db=nan(4,Nntwk,Nntwk,Nfc); Pb=nan(4,Nntwk,Nntwk,Nfc);
  iB=3; PLTLBL=cell(1,4);
  for xx=2
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-'); PLTLBL{xx}=CMstr; % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          label_comm_models=S.delta_stats_all.meta.label_comm_models;  % opts=S.meta.opts;
          deltaR2   = S.delta_deltaR2; 
          stats_all = S.delta_stats_all;
        % Extract data
          Db(xx,:,:,:)=stats_all.(label_level{iB}).(label_targ);  Pb(xx,:,:,:)=stats_all.(label_level{iB}).(label_targ2); % route      
  end
% Check p-vals
  if any([Pb(:)]>=0.05), disp('Check pvals!!'); end
% Plot
  CMAX=0.2; %CMAX = prctile(abs([Db(:)]),97);
% Both
  for xx=1:2
      fig = myfig([PLTLBL{xx} ' (' tmpstr '): \Delta \DeltaR^2'],[-2559 54 563 300+260*4]);
      tl = tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
      for k=1:Nfc
          Dt = squeeze(Db(xx,:,:,k)); ax = nexttile(tl,k);
          imagesc(ax, Dt); axis(ax,'square'); box(ax,'on'); title(ylbl{k});
          add_mat_box(ax, Dt, 'Color','k','LineWidth',0.5,'Outer',true);
          clim(ax, [-CMAX CMAX]); colormap(ax,cm2); cb=colorbar(ax); font(fntsz-4,usefont);
          set(ax,'YTick',1:Nntwk,'YTickLabels',LBL_ntwk,'XTick',1:Nntwk,'XTickLabels',LBL_ntwk); xtickangle(40); ytickangle(0);
      end
  end

% Save
  svstr='2_ntwk_';
  set(gcf,'Renderer','painters'); exportgraphics(gcf,['~/Downloads/' svstr tmptag '3-both.pdf'],'ContentType','vector');



% --- Level 1: Route & Diff separate
  % label_targ='ntwk_deltaR2'; label_targ2='ntwk_p'; tmpstr=[use_myelin_str ' - Caliber'];
  for xx=[1 4]
          CC=testset{xx}; CMstr=strjoin(CMLBL_simple(CC),'-');                               % data
        % Load
          svdir=[locDeriv '/x_dual_BICS/grp_' parc]; load([svdir '/' tmptag CMstr '.mat']);
          label_comm_models=S.delta_stats_all.meta.label_comm_models;  % opts=S.meta.opts;
          deltaR2   = S.delta_deltaR2; 
          stats_all = S.delta_stats_all;

        % Plot
          for k=1:2
              opts.viz.usek=k;
              figs = bicsdual.viz_BICS_dual_allFC_summary(deltaR2, stats_all, ylbl, pinfo, opts.viz);
              ylabel(sprintf('\\mu %s%s','\Delta', '\DeltaR^2'));
          end
  end






% SCRATCH BELOW THIS LINE




% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % %% =======   Modeling 1b: Hiearchical: con & dis separate
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % -------------------------------------------------------------------------
% % % % + Hiearchical regression vs binary base model
% % % % + 2nd-level confirmation by residualization
% % % % + Looped over CON & DIS
% % % 
% % % 
% % % % -------------------------------------------------------------------------
% % % %% ==========================      Global      ============================
% % % % -------------------------------------------------------------------------
% % % 
% % %   Nit=2; Ntest=2;                    % loop num
% % %   CC=5;                              % Switch for communication model
% % %   useCMLBL=[CMLBL_simple{CC}];
% % % 
% % % % Euclidean distance
% % %   Xd=EDnrm; DISTLBL='ED'; Xdz=nzzscore(Xd,0);
% % % 
% % % % --- Communication models
% % % % % % % caliber & myelin
% % % % % %   X1z=ALLCM{1}{CC}; X1_str=CMLBLnx{1}{CC}; X1_str=strrep(X1_str,'-','');    % caliber
% % % % % %   X2z=ALLCM{2}{CC}; X2_str=CMLBLnx{2}{CC}; X2_str=strrep(X2_str,'-','');    % myelin
% % % % myelin & delay
% % %   X1z=ALLCM{2}{CC};    X1_str=CMLBLnx{2}{CC};    X1_str=strrep(X1_str,'-','');    % myelin
% % %   X2z=ALLCM_delay{CC}; X2_str=CMLBLnx_delay{CC}; X2_str=strrep(X2_str,'-','');    % delay
% % % 
% % % 
% % % % Communication model for binary
% % %   Xbz=ALLCM_bin{CC}; Xb_str=CMLBLnx_bin{CC}; Xb_str=strrep(Xb_str,'-','');    % binary
% % % 
% % % % TMP PLOTS
% % %   % % % CLIM=[-3 3];
% % %   % % % connplot(X1z,X1_str,pinfo); clim(CLIM); colormap(cm2); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk)
% % %   % % % connplot(X2z,X2_str,pinfo); clim(CLIM); colormap(cm2); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk)
% % %   % % % connplot(Xbz,Xb_str,pinfo); clim(CLIM); colormap(cm2); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk)
% % %   % % % connplot(Xdz,DISTLBL,pinfo); clim(CLIM); colormap(cm2); set(gca,'XTickLabels',LBL_ntwk,'YTickLabels',LBL_ntwk)
% % % 
% % % 
% % % % Global model variables
% % %   deltaR2_g_ALL=zeros(Nit*Ntest,Nfc);
% % %   pval_g_ALL=zeros(Nit*Ntest,Nfc);
% % %   R2_g_rsd_1_ALL=zeros(Nit*Ntest,Nfc);
% % %   R2_g_rsd_2_ALL=zeros(Nit*Ntest,Nfc);
% % %   R2_g_rsd_3_ALL=zeros(Nit*Ntest,Nfc);
% % %   % % % R2_g_ALL=cell(1,Nit);   
% % %   % % % corr_pred_g_ALL=cell(1,Nit);
% % % 
% % % % Network model variables
% % %   deltaR2_ntwk_ALL=zeros(Nit*Ntest,Nntwk,Nntwk,Nfc);
% % %   pval_ntwk_ALL=zeros(Nit*Ntest,Nntwk,Nntwk,Nfc);
% % %   R2_ntwk_rsd_1_ALL=zeros(Nit*Ntest,Nntwk,Nntwk,Nfc);
% % %   R2_ntwk_rsd_2_ALL=zeros(Nit*Ntest,Nntwk,Nntwk,Nfc);
% % %   R2_ntwk_rsd_3_ALL=zeros(Nit*Ntest,Nntwk,Nntwk,Nfc);
% % % 
% % % % Metadata
% % %   PLTLBL_ALL=cell(Nit*Ntest,1);
% % %   Dcdastr_ALL=cell(1,Nit);
% % % 
% % %   it=0;
% % % 
% % % %% LOOP OVER CON & DIS
% % % for ss=1:Nit
% % % 
% % %     % Combine all X predictors
% % %       useD_base={Xbz      Xdz};         useD_test={X1z      X2z};
% % %       DLBL_base={Xb_str   DISTLBL};     DLBL_test={X1_str   X2_str};
% % %       useN_b=length(useD_base);         useN_t=length(useD_test); 
% % %       ixd=useN_b + 1; 
% % %       iy=ixd+1; 
% % % 
% % %     % Filter for ALL vs CON vs DIS edges
% % %       switch ss 
% % %           case 1 
% % %               Dcdastr=' (CON)'; 
% % %               for kk=1:useN_b; useD_base{kk}(~masknz)=0; end
% % %               for kk=1:useN_t; useD_test{kk}(~masknz)=0; end
% % %           case 2 
% % %               Dcdastr=' (DIS)'; 
% % %               for kk=1:useN_b; useD_base{kk}(~maskdis)=0; end
% % %               for kk=1:useN_t; useD_test{kk}(~maskdis)=0; end
% % %       end    
% % % 
% % % 
% % %   % ----- LOOP OVER TEST PREDICTORS
% % % 
% % %     for tt = 1 : useN_t
% % % 
% % %         it=it+1;
% % % 
% % %         % ---- Combine & compute metadata
% % %         % Combine X & Y
% % %           ALLD=cell(1,iy);            ALLD(1:useN_b)=useD_base;
% % %           ALLD{ixd}=useD_test{tt};    ALLD{iy}=Yz_m; 
% % %           PLTLBL_ALL{it}=[DLBL_test{tt} Dcdastr];
% % % 
% % %         % --- Variables for this loop
% % %         % Global model variables
% % %           pval_g=zeros(1,Nfc);
% % %           deltaR2_g=zeros(1,Nfc);
% % %           R2_g_rsd_1=zeros(1,Nfc);
% % %           R2_g_rsd_2=zeros(1,Nfc);
% % %           R2_g_rsd_3=zeros(1,Nfc);
% % %           % % % R2_g=zeros(2,Nfc);   
% % %           % % % corr_pred_g=zeros(2,Nfc);
% % % 
% % %         % Network model variables
% % %           pval_ntwk=zeros(Nntwk,Nntwk,Nfc);
% % %           deltaR2_ntwk=zeros(Nntwk,Nntwk,Nfc);
% % %           R2_ntwk_rsd_1=zeros(Nntwk,Nntwk,Nfc);
% % %           R2_ntwk_rsd_2=zeros(Nntwk,Nntwk,Nfc);
% % %           R2_ntwk_rsd_3=zeros(Nntwk,Nntwk,Nfc);
% % % 
% % %         disp([' --> Running test for: ' PLTLBL_ALL{it}])
% % % 
% % %         % ---- Global modeling
% % %           for ii = 1 : Nfc
% % %             % Extract all X data with the FC data for this round
% % %               D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
% % %             % Remove 0-valued edges
% % %               Dt=rmzeros(D); X=Dt(:,1:ixd); y=Dt(:,iy);
% % %             % Modeling: BASE
% % %               X_base=X(:,1:useN_b);
% % %               % % % M=fitlm(X_base,y);
% % %               % % % R2_g(1,ii)=M.Rsquared.Adjusted; 
% % %               % % % corr_pred_g(1,ii)=corr(y, M.Fitted);
% % %             % Modeling: FULL
% % %               % % % M=fitlm(X,y);
% % %               % % % R2_g(2,ii)=M.Rsquared.Adjusted; 
% % %               % % % corr_pred_g(2,ii)=corr(y, M.Fitted);
% % %             % Modeling: Hiearchical
% % %               [deltaR2_g(ii),pval_g(ii)]=compare_models(y,X_base,X);
% % %             % --- Orthogonalization tests
% % %               X_test=X(:,useN_b+1);
% % %             % 1. Regress binary & ED out of test predictor & rerun model
% % %               M=fitlm(X_base,X_test);
% % %               X_resid=M.Residuals.Raw;
% % %               M=fitlm(X_resid,y);
% % %               R2_g_rsd_1(ii)=M.Rsquared.Adjusted;
% % %             % 2. Regress binary out of test predictor & rerun model with ED
% % %               M=fitlm(X_base(:,1),X_test);
% % %               X_resid=M.Residuals.Raw;
% % %               M=fitlm([X_resid X_base(:,2)],y);
% % %               R2_g_rsd_2(ii)=M.Rsquared.Adjusted;
% % %             % 3. Regress binary out of test predictor & rerun model without ED
% % %               M=fitlm(X_base(:,1),X_test);
% % %               X_resid=M.Residuals.Raw;
% % %               M=fitlm(X_resid,y);
% % %               R2_g_rsd_3(ii)=M.Rsquared.Adjusted;
% % %           end
% % %         % Store outputs
% % %           deltaR2_g_ALL(it,:)=deltaR2_g; 
% % %           pval_g_ALL(it,:)=pval_g;
% % %           R2_g_rsd_1_ALL(it,:)=R2_g_rsd_1;
% % %           R2_g_rsd_2_ALL(it,:)=R2_g_rsd_2;
% % %           R2_g_rsd_3_ALL(it,:)=R2_g_rsd_3;
% % % 
% % % 
% % % 
% % %         % ========== Pairwise networks  
% % % 
% % %           for ii = 1 : Nfc
% % %             % Extract all X data with the FC data for this round
% % %               D=cell(1,iy);   D(1:ixd)=ALLD(1:ixd); D{iy}=ALLD{iy}{ii};
% % %             % Specify model
% % %               disp(['Running model for y-data: ' ylbl{ii}])
% % % 
% % %               for iio = 1 : Nntwk
% % %                   nio=pinfo.cis{iio};
% % %                   for iii = iio : Nntwk
% % %                       nii=pinfo.cis{iii};
% % %                       ALLDtmp=cell(size(D)); for ii2=1:iy; ALLDtmp{ii2}=D{ii2}(nio,nii); end 
% % %                       Dt=rmzeros(ALLDtmp); X=Dt(:,1:ixd); y=Dt(:,iy);
% % % 
% % % 
% % %                       % Modeling: BASE
% % %                       X_base=X(:,1:useN_b);
% % %                       % % % M=fitlm(X_base,y);
% % %                       % % % R2_g(1,ii)=M.Rsquared.Adjusted; 
% % %                       % % % corr_pred_g(1,ii)=corr(y, M.Fitted);
% % %                     % Modeling: FULL
% % %                       % % % M=fitlm(X,y);
% % %                       % % % R2_g(2,ii)=M.Rsquared.Adjusted; 
% % %                       % % % corr_pred_g(2,ii)=corr(y, M.Fitted);
% % %                     % Modeling: Hiearchical
% % %                       [deltaR2_ntwk(iio,iii,ii),pval_ntwk(iio,iii,ii)]=compare_models(y,X_base,X);
% % %                     % % % % --- Orthogonalization tests
% % %                     % % %   X_test=X(:,useN_b+1);
% % %                     % % % % 1. Regress binary & ED out of test predictor & rerun model
% % %                     % % %   M=fitlm(X_base,X_test);
% % %                     % % %   X_resid=M.Residuals.Raw;
% % %                     % % %   M=fitlm(X_resid,y);
% % %                     % % %   R2_ntwk_rsd_1(iio,iii,ii)=M.Rsquared.Adjusted;
% % %                     % % % % 2. Regress binary out of test predictor & rerun model with ED
% % %                     % % %   M=fitlm(X_base(:,1),X_test);
% % %                     % % %   X_resid=M.Residuals.Raw;
% % %                     % % %   M=fitlm([X_resid X_base(:,2)],y);
% % %                     % % %   R2_ntwk_rsd_2(iio,iii,ii)=M.Rsquared.Adjusted;
% % %                     % % % % 3. Regress binary out of test predictor & rerun model without ED
% % %                     % % %   M=fitlm(X_base(:,1),X_test);
% % %                     % % %   X_resid=M.Residuals.Raw;
% % %                     % % %   M=fitlm(X_resid,y);
% % %                     % % %   R2_ntwk_rsd_3(iio,iii,ii)=M.Rsquared.Adjusted;
% % % 
% % %                   end
% % %               end
% % %           end
% % % 
% % %         % Store outputs
% % %           deltaR2_ntwk_ALL(it,:,:,:)=deltaR2_ntwk; 
% % %           pval_ntwk_ALL(it,:,:,:)=pval_ntwk;
% % %           % % % R2_ntwk_rsd_1_ALL(it,:,:,:)=R2_ntwk_rsd_1;
% % %           % % % R2_ntwk_rsd_2_ALL(it,:,:,:)=R2_ntwk_rsd_2;
% % %           % % % R2_ntwk_rsd_3_ALL(it,:,:,:)=R2_ntwk_rsd_3;
% % %     end
% % %       Dcdastr_ALL{ss}=Dcdastr;
% % % end
% % % 
% % % 
% % % 
% % % %% ===== Visualization
% % % 
% % % % =======  GLOBAL: myelin & caliber colored
% % % 
% % %       str2=['Global model, ' useCMLBL str];  
% % %       % cmap_tmp=[cmap_pred(1:2,:); cmap_pred(1:2,:)]; % caliber myelin
% % %       % cmap_tmp=[cmap_pred(1,:); cmap_intract(1,:); cmap_pred(1,:); cmap_intract(1,:)]; % caliber delay
% % %       cmap_tmp=[cmap_pred(2,:); cmap_intract(1,:); cmap_pred(2,:); cmap_intract(1,:)]; % myelin delay
% % %     % -- deltaR2
% % %       ppostmp=[-2559 541 1132 256];
% % %       D1=deltaR2_g_ALL; D2=pval_g_ALL; D1(D2>0.05)=0; YLIM=ceil(max(D1(:))*100)/100;
% % %       myfig(['1+' strjoin(DLBL_base,'+') '+X: ' str2],ppostmp); B=bar(D1',0.7); font(fntsz,usefont);
% % %       ylabel('\DeltaR^2','FontWeight','bold'); % title(['Base: FC ~ 1 + ' strjoin(DLBL_base,' + ')]); 
% % %       ii=1; B(ii).FaceColor=cmap_tmp(ii,:);  B(ii).FaceAlpha=1;
% % %       ii=2; B(ii).FaceColor=cmap_tmp(ii,:);  B(ii).FaceAlpha=1;
% % %       ii=3; B(ii).FaceColor=cmap_tmp(ii,:);  B(ii).FaceAlpha=.5;
% % %       ii=4; B(ii).FaceColor=cmap_tmp(ii,:);  B(ii).FaceAlpha=.5;
% % %       grid on; set(gca,'XGrid','off'); ylim([0 YLIM]);  
% % %       set(gca,'XTickLabels',ylbl,'XTickLabelRotation',20);
% % %       legend(PLTLBL_ALL,'Location','bestoutside');
% % % 
% % % 
% % % 
% % % % =======  GLOBAL: gray
% % % 
% % %       str2=['Global model: vs binary, ' useCMLBL str];  cmap_tmp=gray(it);
% % %     % -- deltaR2
% % %       ppostmp=[-2559 535 1369 262];
% % %       D1=deltaR2_g_ALL; D2=pval_g_ALL; D1(D2>0.05)=0; YLIM=ceil(max(D1(:))*100)/100;
% % %       myfig([str2 str],ppostmp); B=bar(D1',0.5); font(fntsz,usefont);
% % %       ylabel('\DeltaR^2','FontWeight','bold'); title(['Base: FC ~ 1 + ' strjoin(DLBL_base,' + ')]); 
% % %       for bb=1:it; B(bb).FaceColor=cmap_tmp(bb,:); end 
% % %       grid on; set(gca,'XGrid','off'); ylim([0 YLIM]);  
% % %       set(gca,'XTickLabels',ylbl,'XTickLabelRotation',0);
% % %       legend(PLTLBL_ALL,'Location','bestoutside');
% % % 
% % %     % -- Residualization 1: binary & ED regressed out
% % %       ppostmp=[-2559 535 1369 262];
% % %       D1=R2_g_rsd_1_ALL; D2=pval_g_ALL; D1(D2>0.05)=0; YLIM=ceil(max(D1(:))*100)/100;
% % %       myfig([str2 str],ppostmp); B=bar(D1',0.5); font(fntsz,usefont);
% % %       ylabel('\DeltaR^2','FontWeight','bold'); title(['Base: FC ~ 1 + ' strjoin(DLBL_base,' + ')]); 
% % %       for bb=1:it; B(bb).FaceColor=cmap_tmp(bb,:); end 
% % %       grid on; set(gca,'XGrid','off'); ylim([0 YLIM]);  
% % %       set(gca,'XTickLabels',ylbl,'XTickLabelRotation',0);
% % %       legend(PLTLBL_ALL,'Location','bestoutside');
% % % 
% % % 
% % % 
% % %  % =======  NETWORK LEVEL (CON & DIS separate)
% % % 
% % % 
% % %       useppos=[498 170 -155 -481];    
% % %     % ----------- Matrix Plots
% % %       axopt={'XTick',1:Nntwk,'XTickLabel',LBL_ntwk,'XTickLabelRotation',40, ...
% % %              'YTick',1:Nntwk,'YTickLabel',LBL_ntwk,'YTickLabelRotation',40,'FontSize',fntsz};
% % %       meshopt={'EdgeColor','k','LineWidth',1.5};
% % % 
% % %       for it=1:4
% % %           str2=['NTWK' PLTLBL_ALL{it}]; ppostmp2=getPlotPos([-2559 useppos(it) 315 299],Nfc,0);
% % %         % --- delta R2
% % %           D1=deltaR2_ntwk_ALL; D2=pval_ntwk_ALL; D1(D2>0.05)=0;
% % %           D=squeeze(D1(it,:,:,:));  TMPSTR='\DeltaR^2: '; usecmap=cm3; CLIM=[0 ceil(max(D(:))*10)/10];
% % %           TMPSTR2=[TMPSTR ', ' str2]; 
% % %           for ii=1:Nfc
% % %             % Matrix with blank upper tri
% % %               D2=D(:,:,ii)'; tmask=tril(true(size(D2))); myfig([TMPSTR2 ', ' ylbl{ii}],ppostmp2(ii,:));
% % %               imagesc(D2,'AlphaData',tmask); colormap(usecmap); colorbar; axis square; clim(CLIM); font(fntsz,usefont);
% % %               set(gca,axopt{:}); font(fntsz,usefont); hold on; title(ylbl{ii}); 
% % %               for I=1:Nntwk; for J=1:I; rectangle('Position',[J-0.5,I-0.5,1,1],meshopt{:}); end; end % draw box around all data in lower tri + main diag
% % %               axis off;
% % %               % set(gca,'XTickLabel','','YTickLabel',''); title ''; % colorbar off;
% % %           end
% % %       end
% % % 
% % % 
% % % 
% % % % % % % =======  NETWORK LEVEL (CON & DIS together)
% % % % % % 
% % % % % % % Setup
% % % % % %   str2='NTWK: ';
% % % % % %   MDLSPECGEN=strjoin(PLTLBL_hi, '+');
% % % % % %   TMPSTR2=[useCMLBL ' w/ xtion'];
% % % % % %   TMPLBL3={'SC>0' 'SC=0'};
% % % % % %   [I,J]=find(triu(true(Nntwk)));          % row & col indices for upper tri + main diag
% % % % % %   idx=sub2ind([Nntwk, Nntwk],I,J);        % linear indices
% % % % % %   Nntwkpairs=numel(idx);
% % % % % %   boxopt1={'FactorSeparator',1,'Colors','k','Symbol','','Widths',0.6};
% % % % % %   boxopt2={'filled','MarkerFaceAlpha',0.6,'MarkerEdgeColor','k','LineWidth',0.3};
% % % % % %   boxopt3={'FaceAlpha',0.25,'LineWidth',1.5};
% % % % % %   ppostmp=[-2559 445 1357 352];
% % % % % % 
% % % % % % % ----- R2 values: box plots
% % % % % % % Reshape R2 values (network-pairs, Nit, FC)
% % % % % %   D=R2_ntwk_ALL; TMPSTR='R^2';
% % % % % %   DPLT=zeros(Nntwkpairs,Nit,Nfc);
% % % % % %   for xx=1:Nfc
% % % % % %       for ss=1:Nit
% % % % % %         tmp=D{ss}(:,:,xx); DPLT(:,ss,xx)=tmp(idx);
% % % % % %       end
% % % % % %   end
% % % % % %   DPLTvec    = reshape(DPLT,[],1);                                        % reshape for boxplot [29*6*2 x 1]
% % % % % %   iDgrps_fc  = repelem(1:Nfc, Nntwkpairs*Nit)';                           % Indices for FC groupings (outer group)
% % % % % %   iDgrps_cm  = repmat(repelem(1:Nit, Nntwkpairs), [1 Nfc])';              % Indices for CON-DIS groupings (inner group)
% % % % % %   iDgrps_cm  = iDgrps_cm(:);
% % % % % % 
% % % % % % % -- Boxplot
% % % % % %   myfig([str2 TMPSTR2 str],ppostmp); hold on;
% % % % % %   boxplot(DPLTvec,{iDgrps_fc,iDgrps_cm},boxopt1{:});
% % % % % %   set(gca,'XTickLabel',''); ylabel(TMPSTR); font(fntsz,usefont);
% % % % % %   grid on; set(gca,'XGrid','off','XMinorGrid','off');
% % % % % % % -- Overlay data points
% % % % % %   h=findobj(gca,'Tag','Box'); Nboxes=Nit*Nfc; tpos=nan(Nboxes,1);          % compute x positions for each box
% % % % % %   for i=1:Nboxes; boxX=get(h(i),'XData'); tpos(Nboxes-i+1)=mean(boxX); end % reverse order
% % % % % % % Plot each set of points
% % % % % %   ptIdx = 0;
% % % % % %   for xx = 1:Nfc
% % % % % %       for ss = 1:Nit
% % % % % %           ptIdx=ptIdx+1;
% % % % % %           x=repmat(tpos(ptIdx),Nntwkpairs,1);
% % % % % %           y=DPLT(:,ss,xx);
% % % % % %           scatter(x+0.05*randn(size(x)),y,50,cmap_discon(ss,:),boxopt2{:});
% % % % % %       end
% % % % % %   end
% % % % % % % Set tick labels
% % % % % %   xticks(mean(reshape(tpos,Nit,Nfc),1)); xticklabels(ylbl); box on;
% % % % % % % Mod boxes
% % % % % % for i = 1:length(h)
% % % % % %       ibox=h(i);                                                         % this box
% % % % % %       icm=mod(length(h)-i,Nit)+1; boxColor=cmap_discon(icm, :);             % Boxes are stored in reverse order
% % % % % %       x=get(ibox,'XData'); y=get(ibox,'YData');                          % Get X and Y data for patching
% % % % % %       patch(x,y,boxColor,'EdgeColor',boxColor,boxopt3{:});                                    % color, thicker lines & lower alpha
% % % % % % end
% % % % % % % Add legend
% % % % % %   hold on; tmp=nan(1,Nit); 
% % % % % %   for ss=1:Nit; tmp(ss)=scatter(nan,nan,50,cmap_discon(ss,:),'filled'); end % blank scatter plot
% % % % % %   legend(tmp,TMPLBL3,'Location','bestoutside');
% % % 



% -------------------------------------------------------------------------









