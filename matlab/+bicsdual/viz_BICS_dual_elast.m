function figs = viz_BICS_dual_elast(elast_all, stats_all, ylbl, pinfo, opts)
% VIZ_BICS_DUAL_ELAST  Elasticity visualizations for dual-CM workflow.
% Handles the M=1 special case (across present levels), and M>1 fallback.
%
% Inputs
%   elast_all : struct with fields like Delta_all (e.g., .L1_route, .L1_diff, .L2_both, .L3_int)
%               Each field is a B x M x Nfc cube of elasticity values.
%   stats_all : struct mirroring elast_all with .global_deltaR2, .model_labels, .block_map, etc.
%   ylbl      : 1xNfc band labels
%   pinfo     : struct with .cis and .clabels_short (for axes)
%   opts      : struct (all optional)
%       .norm             : 'none' (default) | 'by_delta'     % elementwise (/ max(Delta, floor))
%       .norm_floor       : 1e-6 (default)
%       .fc_agg           : 'mean' (default) | 'median'
%       .figpos_meq1      : [x y w h] for special-case heatmaps
%       .ribbons_units    : auto from .elasticity_mode, else char
%       .elasticity_mode  : 'zadd' (default) | 'xscale'  % only for axis label text
%       .ribbons_xdim     : 'bands' (default here) | 'models'
%
% 2025 Mark C Nelson (MNI)
% -------------------------------------------------------------------------

if nargin<5, opts = struct; end
normMode   = getopt(opts,'norm','none');                    % 'none' | 'by_delta'
normFloor  = getopt(opts,'norm_floor',1e-6);
fc_agg     = getopt(opts,'fc_agg','mean');                  % 'mean'|'median'
fnt        = 16;
figpos1    = getopt(opts,'figpos1',[80 80 900 750]);
ribbons_x  = getopt(opts,'ribbons_xdim','bands');           % default 'bands' for elasticity
emode      = getopt(opts,'elasticity_mode','zrow');         % just for y-axis label text
labelON    = getopt(opts,'labelON',true);
lineopt    = {'LineWidth',2,'MarkerSize',10};
labels_CMs = stats_all.meta.label_comm_models;
labels_int = stats_all.meta.interaction_label;

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

% Gather levels present + M per level
fields = {'L1_route','L1_diff','L2_both','L3_int'};
nice   = {'L1-Route','L1-Diff','L2-Both','L3-Both+Int'};

present = false(1,numel(fields));
M_per   = nan(1,numel(fields));
for k=1:numel(fields)
    f = fields{k};
    if isfield(elast_all,f) && ~isempty(elast_all.(f))
        present(k) = true;
        [~,M_per(k),~] = size(elast_all.(f));
    end
end

% Optional norm by Delta (elementwise)
if strcmpi(normMode,'by_delta')
    assert(isfield(opts,'Delta_all') && ~isempty(opts.Delta_all), ...
        'norm="by_delta" requires opts.Delta_all with matching fields.');
    Delta_all = opts.Delta_all;
else
    Delta_all = struct(); % dummy
end

YLBL=elasticity_ylabel(normMode, emode);

% ============ SPECIAL CASE: all present have M==1 ============
if any(present) && all(M_per(present)==1)
    % Compute FC-aggregated block maps and common clims
    Dmats = {}; Emax = 0; lev_names = {};
    Smats = {};                                         % kept in case want to add SI-like thing
    nlev  = sum(present);

    for k=find(present)
        f  = fields{k};
        St = stats_all.(f);
        blk_ij = St.block_map; Nntwk = St.Nntwk;

        E = elast_all.(f); % B x 1 x Nfc
        if strcmpi(normMode,'by_delta')
            D = opts.Delta_all.(f);
            denom = max(D, normFloor);
            E = E ./ denom;
        end
        Eagg = aggfun(E);    % B x 1

        Em = blocks_to_mat(Eagg(:,1), blk_ij, Nntwk, true);
        Dmats{end+1} = Em; %#ok<AGROW>
        % Emax = max(Emax, prctile( abs(Em(:)),95 ));
        Emax = max(Emax, prctile( abs(Em(:)),95));
        lev_names{end+1} = nice{k}; %#ok<AGROW>
    end
    % if Emax<=0, Emax = 1; end

    % ---- Heatmaps (rows = levels, common clims)
    figs.elast_heatmaps_Meq1 = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect ': ' YLBL],figpos1);
    tl = tiledlayout(nlev,1,'TileSpacing','compact','Padding','compact');
    for r = 1:nlev
        f  = fields{find(present, r, 'first')}; % just to fetch Nntwk each row
        Nntwk = stats_all.(f).Nntwk;
        nexttile(tl,r);
        imagesc(Dmats{r}); axis square; box on; colorbar;
        clim([-Emax Emax]); tick_rsn(pinfo,Nntwk); set(gca,'FontSize',fnt+2);
        if labelON
            title(sprintf('%s — Elasticity', lev_names{r}));
        end
    end
    % tl.Title.String   = elasticity_title(normMode, emode);
    % tl.Title.FontSize = fnt+2;
    colormap(smartcmaps('bentcoolwarm'));

    % ---- Ribbons overlay (lines = levels, x = bands)
    figs.elast_ribbons_Meq1 = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect ': ' YLBL],[120 120 1000 380]); hold on;
    cols = lines(nlev);
    x = 1:numel(ylbl);
    ridx = 0;
    for k=find(present)
        ridx = ridx+1;
        f = fields{k};
        E = elast_all.(f);   % B x 1 x Nfc
        if strcmpi(normMode,'by_delta')
            D = opts.Delta_all.(f);
            denom = max(D, normFloor);
            E = E ./ denom;
        end
        % Median over blocks per band
        Eb = squeeze(median(E,1,'omitnan'));  % 1 x 1 x Nfc -> 1 x Nfc
        plot(x,Eb,':o','Color',cols(ridx,:),'MarkerFaceColor',cols(ridx,:),lineopt{:});
    end
    grid on; box on;
    xlim([0.7 numel(ylbl)+0.3]); xticks(x); xticklabels(ylbl); xtickangle(0);
    % ylabel(elasticity_ylabel(normMode, emode));
    xlabel('FC bands'); legend(lev_names, 'Location','northeastoutside');
    title('Elasticity across bands'); set(gca,'FontSize',fnt+2);
    return;
end

% ============ FALLBACK: per-level small multiples (M>1 in any level) ====
figs = struct();

% Per-level grids and ribbons
for k=1:numel(fields)
    f = fields{k};
    if ~isfield(elast_all,f), continue; end
    E = elast_all.(f); % B x M x Nfc
    St = stats_all.(f);
    blk_ij = St.block_map; Nntwk = St.Nntwk;
    [B,M,Nfc] = size(E); %#ok<ASGLU>

    % Optional normalization
    if strcmpi(normMode,'by_delta')
        D = opts.Delta_all.(f);
        denom = max(D, normFloor);
        E = E ./ denom;
    end

    % Heatmaps: 1 row (elasticity), M columns, FC-aggregated
    Eagg = aggfun(E);                  % B x M
    emax = max(abs(Eagg(:))); if emax<=0, emax=1; end

    figs.(f) = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect ': ' YLBL],getopt(opts,'figpos',[100 100 320+260*M 360]));
    tl = tiledlayout(1,M,'TileSpacing','compact','Padding','compact');
    for m=1:M
        nexttile(tl,m);
        Em = blocks_to_mat(Eagg(:,m), blk_ij, Nntwk, true);
        imagesc(Em); axis square; box on; colorbar;
        title(sprintf('Elasticity | %s', St.model_labels{m}));
        clim([-emax emax]); tick_rsn(pinfo,Nntwk); set(gca,'FontSize',fnt+2);
    end
    tl.Title.String = [nice{k} ' — Elasticity']; tl.Title.FontSize = fnt+2;

    % Ribbons: bands on x, lines = models
    figs.([f '_ribbons']) = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect ': ' YLBL],[120 120 1000 380]);
    G = squeeze(median(E,1,'omitnan'));    % M x Nfc
    plot(G.',':o',lineopt{:}); grid on; box on;
    xticks(1:size(G,2)); xticklabels(ylbl); xtickangle(0);
    legend(St.model_labels,'Location','northeastoutside');
    % xlabel('FC bands'); ylabel(elasticity_ylabel(normMode, emode));
    title([nice{k} ' — Elasticity']); set(gca,'FontSize',fnt+2);
end

% -------------------------------------------------------------------------
end

% ----------------- helpers -----------------
function v = getopt(s, name, default)
if ~isstruct(s), v = default; return; end
if isfield(s,name), v=s.(name); else, v=default; end
end

function t = elasticity_title(normMode, emode)
switch lower(normMode)
    case 'none'
        switch lower(emode)
            case 'zadd',  t = 'Elasticity (ΔR^2 per SD; FC-agg)';
            case 'xscale',t = 'Elasticity (per log scale; FC-agg)';
            otherwise,    t = 'Elasticity (FC-agg)';
        end
    case 'by_delta'
        switch lower(emode)
            case 'zadd',  t = 'Normalized Elasticity (per SD / ΔR^2; FC-agg)';
            case 'xscale',t = 'Normalized Elasticity (per log scale / ΔR^2; FC-agg)';
            otherwise,    t = 'Normalized Elasticity (/ΔR^2; FC-agg)';
        end
    otherwise
        t = 'Elasticity';
end
end

function yl = elasticity_ylabel(normMode, emode)
switch lower(normMode)
    case 'none'
        switch lower(emode)
            case 'zrow',  yl = 'Elasticity (\DeltaR^2 per SD, row-wise)';
            case 'zadd',  yl = 'Elasticity (\DeltaR^2 per SD)';
            case 'xscale',yl = 'Elasticity (\DeltaR^2 per log scale)';
            otherwise,    yl = 'Elasticity';
        end
    case 'by_delta'
        switch lower(emode)
            case 'zrow',  yl = 'Norm. Elasticity (% per SD, row-wise)';
            case 'zadd',  yl = 'Norm. Elasticity (% per SD)'; 
            case 'xscale',yl = 'Norm. Elasticity (% per log scale)';
            otherwise,    yl = 'Norm. Elasticity';
        end
    otherwise
        yl = 'Elasticity';
end
end

function M = blocks_to_mat(bvec, blk_ij, Nntwk, symmetrize)
M = zeros(Nntwk,Nntwk);
for b=1:size(blk_ij,1)
    i=blk_ij(b,1); j=blk_ij(b,2);
    M(i,j) = bvec(b);
    if symmetrize, M(j,i) = bvec(b); end
end
end

function tick_rsn(pinfo,Nntwk)
xticks(1:Nntwk); yticks(1:Nntwk);
if isfield(pinfo,'clabels_short') && numel(pinfo.clabels_short)==Nntwk
    xticklabels(pinfo.clabels_short); yticklabels(pinfo.clabels_short);
    xtickangle(40); ytickangle(0);
end
end
