function figs = viz_BICS_dual(Delta_all, SI_all, stats_all, ylbl, pinfo, opts)
% VISUALIZATION (with special handling for the common M=1 case)
% deltaR2 and SI averaged across FC bands.
%
% 2025 Mark C Nelson (MNI)
% -------------------------------------------------------------------------

if nargin<6, opts=struct; end
ribbons_xdim = getopt(opts,'ribbons_xdim','models');     % 'models' | 'bands'
labelON      = getopt(opts,'labelON',true);
fc_agg       = getopt(opts,'fc_agg','mean');                  % 'mean'|'median'

% figpos1      = getopt(opts,'figpos1',[80 80 900 750]);
fnt = 16;
SIbinBounds = [2 10 30 50]; 
SIbinLabels = {'ns' 'weak' 'moderate' 'strong' 'very strong'};

lineopt={'LineWidth',2,'MarkerSize',10};
labels_CMs      = stats_all.meta.label_comm_models;
labels_int      = stats_all.meta.interaction_label;

if isfield(stats_all.meta,'effect_mode') 
    labels_effect = stats_all.meta.effect_mode;
else
    labels_effect = '';
end

% Choose aggregation across FCs
if strcmpi(fc_agg,'median')
    aggfun = @(X) median(X,3);
else
    aggfun = @(X) mean(X,3);
end

% --- gather available levels in order
fields = {'L1_route','L1_diff','L2_both','L3_int'};
nice   = {'L1-Route','L1-Diff','L2-Both','L3-Both+Int'};

present = false(1,numel(fields));
M_per   = nan(1,numel(fields));
for k=1:numel(fields)
    f = fields{k};
    if isfield(Delta_all,f) && ~isempty(Delta_all.(f))
        present(k) = true;
        [~,M_per(k),~] = size(Delta_all.(f));
    end
end

% ====== SPECIAL CASE: all present levels have M==1 ======
if any(present) && all(M_per(present)==1)
    figs = struct();

    % ---- compute common clims across levels for each metric
    Dmats = {}; Smats = {};
    Dmax  = 0;  Smax  = 0;
    for k=find(present)
        f  = fields{k};
        St = stats_all.(f);
        blk_ij = St.block_map; Nntwk = St.Nntwk;

        % FC-aggregate (mean or median)
        Dagg = aggfun(Delta_all.(f));    % B x 1 x Nfc -> B x 1
        SIagg = aggfun(SI_all.(f));    % B x 1
        % Dagg = mean(Delta_all.(f), 3);   % B x 1 x Nfc -> B x 1 after mean
        % SIagg= mean(SI_all.(f),    3);   % B x 1

        % Map to Nntwk x Nntwk matrices
        Dm = blocks_to_mat(Dagg(:,1), blk_ij, Nntwk, true);
        Sm = blocks_to_mat(SIagg(:,1), blk_ij, Nntwk, true);

        Dmats{end+1} = Dm; %#ok<AGROW>
        Smats{end+1} = Sm; %#ok<AGROW>

        % Dmax = max(Dmax, max(Dm(:)));
        % Smax = max(Smax, max(abs(Sm(:))));
        Dmax = max(Dmax, prctile(Dm(:),95));
        % Smax = max(Smax, prctile(abs(Sm(:)),95));

    end
    if Dmax<=0, Dmax = 1; end
    % if Smax<=0, Smax = 1; end
    Dmax = ceil(Dmax*100)/100;  % tidy upper bound

    % Map raw SI values to categorical bins
    % Sbins = rawSI_to_bins(Smats, SIbinBounds);         % signed bins (default)
    Sbins = rawSI_to_bins(Smats, SIbinBounds, false);  % unsigned bins (1..N+1)
    Smax  = length(SIbinBounds)+1;         


    % ---- Figure A: heatmaps (rows = levels, cols = {ΔR², SI})
    nlev = sum(present);
    figs.heatmaps_Meq1 = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect],getopt(opts,'figpos_meq1',[-2559 54 563 300+260*nlev]));
    tl = tiledlayout(nlev,2,'TileSpacing','compact','Padding','compact');

    lev_idx = find(present);
    for r = 1:nlev
        k  = lev_idx(r);
        f  = fields{k};
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
        clim(axL, [0 Dmax]);
        colormap(axL, smartcmaps('inferno'));                     % per-axes colormap
        cbL = colorbar(axL); %#ok<NASGU>
        tick_rsn(pinfo, Nntwk);
        set(axL,'FontSize',fnt);
        if labelON
            title(axL, sprintf('%s — \\DeltaR^2', nice{k}));
        else
            if r~=nlev, set(axL,'XTickLabel',''); end

        end 

    
      % --- (right) SI
        axR = nexttile(tl,(r-1)*2+2);
        imagesc(axR, Sm); axis(axR,'square'); box(axR,'on');
        add_mat_box(axR, Sm, 'Color','k','LineWidth',0.5,'Outer',true);
      % Plot settings
        clim(axR, [0.5 Smax+0.5]); colormap(axR, slanCM('viridis',Smax));
        %clim(axR, [-Smax Smax]); 
        cbR = colorbar(axR);                                    %#ok<NASGU>
        tick_rsn(pinfo, Nntwk);
        set(axR,'FontSize',fnt); cbR.TickLabels=[];
        if labelON
            title(axR, sprintf('%s — SI', nice{k}));
        else
            if r~=nlev, set(axR,'XTickLabel',''); end
            set(axR,'YTickLabel','');
        end
    end
    % tl.Title.String   = 'ΔR^2 and SI';
    % tl.Title.FontSize = fnt+2;

    % ---- Figure B: ribbons (one plot, lines = levels, x = bands)
    figs.ribbons_Meq1 = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect],[120 120 1000 380]); hold on;
    cols = lines(nlev);
    lev_names = {};
    x = 1:numel(ylbl);
    for r = 1:nlev
        k  = lev_idx(r);
        f  = fields{k};
        G  = stats_all.(f).global_deltaR2;   % M x Nfc, but M==1 => 1 x Nfc
        plot(x, G(1,:),':o','Color',cols(r,:),'MarkerFaceColor',cols(r,:),lineopt{:});
        lev_names{end+1} = nice{k}; %#ok<AGROW>
    end
    grid on; box on;
    xlim([0.7 numel(ylbl)+0.3]); xticks(x); xticklabels(ylbl); xtickangle(0);
    ylabel('\DeltaR^2 (global)'); xlabel('FC bands');
    legend(lev_names, 'Location','northeastoutside');
    title('Global \DeltaR^2 across bands'); set(gca,'FontSize',fnt+2);

    return; % done with special case
end

% ====== FALLBACK: your existing per-level small-multiples and ribbons ======
figs = struct();

% Helper to draw a 2xM per-level grid
function fig = draw_level(level_name, Delta_level, SI_level, stats_level)
    [B,M,Nfc] = size(Delta_level); %#ok<NASGU>
    fig = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect],getopt(opts,'figpos',[100 100 320+260*M 600]));
    tl = tiledlayout(2,M,'TileSpacing','compact','Padding','compact');
    blk_ij = stats_level.block_map; Nntwk = stats_level.Nntwk;
    % Agg across FCs
    Dagg = aggfun(Delta_level);
    SIagg = aggfun(SI_level);
    % Dagg = mean(Delta_level,3);
    % SIagg= mean(SI_level,3);
    dmax = ceil(prctile(Dagg(:),95)*10)/10; 
    smax = prctile(abs(SIagg(:)),95);
    for m=1:M
        % ΔR²
        nexttile(tl,m);
        Dm = blocks_to_mat(Dagg(:,m), blk_ij, Nntwk, true);
        imagesc(Dm); axis square; box on; colorbar; title(sprintf('%s | %s','\DeltaR^2', stats_level.model_labels{m}));
        clim([0 dmax]); tick_rsn(pinfo,Nntwk); set(gca,'FontSize',fnt);
        colormap(smartcmaps('inferno'));
        % SI
        nexttile(tl,M+m);
        Sm = blocks_to_mat(SIagg(:,m), blk_ij, Nntwk, true);
        imagesc(Sm); axis square; box on; colorbar; title(sprintf('%s | %s','SI', stats_level.model_labels{m}));
        clim([-smax smax]); tick_rsn(pinfo,Nntwk); set(gca,'FontSize',fnt);
        colormap(smartcmaps('bentcoolwarm'));
    end
    tl.Title.String = level_name; tl.Title.FontSize = fnt+2;
end

for k=1:numel(fields)
    f = fields{k};
    if isfield(Delta_all,f)
        Dl = Delta_all.(f);
        Sl = SI_all.(f);
        St = stats_all.(f);
        figs.(f) = draw_level(nice{k}, Dl, Sl, St);
    end
end

% Also add simple global ribbons per level (kept flexible)
function fig = draw_ribbons(level_name, stats_level, ribbons_xdim, ylbl, fnt)
    fig = myfig([labels_CMs ', ' labels_int],[120 120 1000 380]);
    G = stats_level.global_deltaR2; % M x Nfc
    % lineopt={'LineWidth',2,'MarkerSize',10};
    switch lower(ribbons_xdim)
        case 'models'
            plot(G,':o',lineopt{:}); grid on; box on;
            xticks(1:size(G,1));
            xticklabels(stats_level.model_labels); xtickangle(35);
            legend(ylbl,'Location','northeastoutside');    % bands in legend
            xlabel('Models'); ylabel('\DeltaR^2 (global)');
        case 'bands'
            plot(G.',':o',lineopt{:}); grid on; box on;
            xticks(1:size(G,2));
            xticklabels(ylbl); xtickangle(0);
            legend(stats_level.model_labels,'Location','northeastoutside');  % models in legend
            xlabel('FC bands'); ylabel('\DeltaR^2 (global)');
        otherwise
            error('Unknown ribbons_xdim: %s (use ''models'' or ''bands'')', ribbons_xdim);
    end
    title([level_name ' — Global \DeltaR^2 by band']); set(gca,'FontSize',fnt);
end

for k=1:numel(fields)
    f = fields{k};
    if isfield(stats_all,f)
        figs.([f '_ribbons']) = draw_ribbons(nice{k}, stats_all.(f), ribbons_xdim, ylbl, fnt);
    end
end
% -------------------------------------------------------------------------
end
