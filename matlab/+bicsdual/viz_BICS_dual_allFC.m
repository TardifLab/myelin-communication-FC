function figs = viz_BICS_dual_allFC(Delta_all, SI_all, stats_all, ylbl, pinfo, opts)
% VISUALIZATION
% viz_BICS_dual without cross-FC aggregation (M=1 case only)
%
% 2025 Mark C Nelson (MNI)
% -------------------------------------------------------------------------

if nargin<6, opts=struct; end
labelON     = getopt(opts,'labelON',true);
usek        = getopt(opts,'usek',1);
matchClim   = getopt(opts,'matchClim',false);

% figpos1      = getopt(opts,'figpos1',[80 80 900 750]);
fnt = 16;
SIbinBounds = [2 10 30 50]; 
Nfc         = length(ylbl);

labels_CMs      = stats_all.meta.label_comm_models;
labels_int      = stats_all.meta.interaction_label;

if isfield(stats_all.meta,'effect_mode') 
    labels_effect = stats_all.meta.effect_mode;
else
    labels_effect = '';
end

fields  = {'L1_route','L1_diff','L2_both','L3_int'};
nice    = {'L1-Route','L1-Diff','L2-Both','L3-Both+Int'};
f       = fields{usek};

% ====== SPECIAL CASE: all present levels have M==1 ======
figs = struct();

% ---- Gather data
Dmats = {}; Smats = {};
Dmax  = 0;  Smax  = 0;
for k=1:Nfc
    
    St = stats_all.(f);
    blk_ij = St.block_map; Nntwk = St.Nntwk;

    % FC
    Dt  = Delta_all.(f)(:,1,k);    % B x 1 x Nfc -> B x 1
    SIt = SI_all.(f)(:,1,k);       % B x 1

    % Map to Nntwk x Nntwk matrices
    Dm = blocks_to_mat(Dt(:,1), blk_ij, Nntwk, true);
    Sm = blocks_to_mat(SIt(:,1), blk_ij, Nntwk, true);

    Dmats{end+1} = Dm; %#ok<AGROW>
    Smats{end+1} = Sm; %#ok<AGROW>

    Dmax = max(Dmax, prctile(Dm(:),95));

end
if Dmax<=0, Dmax = 1; end
Dmax = ceil(Dmax*100)/100;  % tidy upper bound

% Map raw SI values to categorical bins
  Sbins = rawSI_to_bins(Smats, SIbinBounds, false);   % unsigned bins (1..N+1)
  Smax  = length(SIbinBounds)+1; 

% Custom colormaps
  cm_si = slanCM('viridis',Smax); cm_si(1,:) = [0 0 0];                     % SI
  cm_r2 = smartcmaps('inferno');                                            % deltaR2



% ---- Figure A: heatmaps (rows = FC, cols = {ΔR², SI})
nlev = ceil(Nfc/2);
figs.heatmaps1 = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect ', ' nice{usek}],getopt(opts,'figpos_meq1',[-2559 54 563 300+260*nlev]));
tl = tiledlayout(nlev,2,'TileSpacing','compact','Padding','compact');

% lev_idx = find(present);
for r = 1:nlev
    St = stats_all.(f);
    Nntwk = St.Nntwk;

  % Pull prepared mats
    Dm = Dmats{r}; 
    Sm = Sbins{r}; %Sm = Smats{r};

  % --- (left) ΔR²
    axL = nexttile(tl,(r-1)*2+1);
    imagesc(axL, Dm); axis(axL,'square'); box(axL,'on');
    add_mat_box(axL, Dm, 'Color','k','LineWidth',0.5,'Outer',true);
  % Plot settings
    if matchClim
        clim(axL, [0 Dmax]);                                    % match clims across FC
    else
        CLIM=ceil(prctile(Dm(Dm~=0),95)*100)/100;               % FC-specific clims
        clim(axL, [0 CLIM]);
    end
    colormap(axL, cm_r2);                                       % per-axes colormap
    cbL = colorbar(axL); %#ok<NASGU>
    tick_rsn(pinfo, Nntwk);
    set(axL,'FontSize',fnt);
    if labelON
        title(axL, sprintf('%s — \\DeltaR^2', ylbl{r}));
    else
        if r~=nlev, set(axL,'XTickLabel',''); end

    end 


  % --- (right) SI
    axR = nexttile(tl,(r-1)*2+2);
    imagesc(axR, Sm); axis(axR,'square'); box(axR,'on');
    add_mat_box(axR, Sm, 'Color','k','LineWidth',0.5,'Outer',true);
  % Plot settings
    clim(axR, [0.5 Smax+0.5]); colormap(axR, cm_si);
    cbR = colorbar(axR);                                    %#ok<NASGU>
    tick_rsn(pinfo, Nntwk);
    set(axR,'FontSize',fnt); cbR.TickLabels=[];
    if labelON
        title(axR, sprintf('%s — SI', ylbl{r}));
    else
        if r~=nlev, set(axR,'XTickLabel',''); end
        set(axR,'YTickLabel','');
    end
end

% ---- Figure A: heatmaps (rows = FC, cols = {ΔR², SI})
nlev2 = Nfc-nlev;
figs.heatmaps2 = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect ', ' nice{usek}],getopt(opts,'figpos_meq1',[-1995 54 563 300+260*nlev]));
tl = tiledlayout(nlev2,2,'TileSpacing','compact','Padding','compact');

for r = 1:nlev2
    k  = r + nlev;
    St = stats_all.(f);
    Nntwk = St.Nntwk;

  % Pull prepared mats
    Dm = Dmats{k}; 
    Sm = Sbins{k}; %Sm = Smats{r};

  % --- (left) ΔR²
    axL = nexttile(tl,(r-1)*2+1);
    imagesc(axL, Dm); axis(axL,'square'); box(axL,'on');
    add_mat_box(axL, Dm, 'Color','k','LineWidth',0.5,'Outer',true);
  % Plot settings
    if matchClim
        clim(axL, [0 Dmax]);                                    % match clims across FC
    else
        CLIM=ceil(prctile(Dm(Dm~=0),95)*100)/100;               % FC-specific clims
        clim(axL, [0 CLIM]);
    end
    colormap(axL, cm_r2);                                       % per-axes colormap
    cbL = colorbar(axL); %#ok<NASGU>
    tick_rsn(pinfo, Nntwk);
    set(axL,'FontSize',fnt);
    if labelON
        title(axL, sprintf('%s — \\DeltaR^2', ylbl{k}));
    else
        if r~=nlev2, set(axL,'XTickLabel',''); end

    end 


  % --- (right) SI
    axR = nexttile(tl,(r-1)*2+2);
    imagesc(axR, Sm); axis(axR,'square'); box(axR,'on');
    add_mat_box(axR, Sm, 'Color','k','LineWidth',0.5,'Outer',true);
  % Plot settings
    clim(axR, [0.5 Smax+0.5]); colormap(axR, cm_si);
    cbR = colorbar(axR);                                    %#ok<NASGU>
    tick_rsn(pinfo, Nntwk);
    set(axR,'FontSize',fnt); cbR.TickLabels=[];
    if labelON
        title(axR, sprintf('%s — SI', ylbl{k}));
    else
        if r~=nlev2, set(axR,'XTickLabel',''); end
        set(axR,'YTickLabel','');
    end
end

% -------------------------------------------------------------------------
end
