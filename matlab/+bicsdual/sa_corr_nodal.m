function [out, allsurfh] = sa_corr_nodal(stats_all, sa_axis, ylbl, pinfo, opts)
% bicsdual.sa_corr_nodal
% Correlate nodal ΔR² maps (from Stage A stats) with a supplied S-A axis.
%
% Inputs
%   stats_all : struct from Stage A with fields like .L1_route, .L1_diff, .L2_both
%               Each has .node_deltaR2 (M x Nnode x Nfc) and .model_labels, etc.
%   sa_axis   : Nnode x 1 vector (S-A axis per node)
%   ylbl      : 1 x Nfc cellstr (band labels, e.g., {'BOLDin','theta',...})
%   pinfo     : struct with .cis, .clabels_short (used only for surfaces)
%   opts      : (all optional)
%       .levels        : cellstr subset of {'L1_route','L1_diff','L2_both'} (default: present)
%       .corr_type     : 'Spearman' (default) | 'Pearson'
%       .absval        : false (if true, use |r|)
%       .figpos_bars   : [x y w h], default [100 100 1100 420]
%       .m_eq1_colors  : [] (use lines) or nlev x 3 RGB double in [0,1]
%       .bands         : indices of bands to include (default: all)
%       .labelON       : true (if false, some plot labels turned off)
%       .make_surfaces : false (default) — if true, draws ΔR² surfaces
%       .surf_bands    : indices for surfaces (default: [])
%       .surf_level    : which level to surface when M=1 case (default: 'L2_both')
%       .surf_title    : char, default 'Cortical Features'
%       .surf_cmax     : [] auto (95th pct of |ΔR²|); or numeric >0
%       .surf_cmap     : [] (use your env default); or an Mx3 colormap
%       .matchClim     : 0/1 match clims across surfaces (default = 1)
%
% Output
%   out.r       : struct with r vectors per level (1 x Nfc when M=1; else M x Nfc)
%   out.p       : same as r with p-values
%   out.fig_bars: bars figure handle
%   out.fig_surf: array of surface fig handles (if any)
%   out.meta    : bookkeeping (levels, bands, corr_type, etc.)
%   allsurfh    : handles to all surf objects
%
% 2025 Mark C Nelson (MNI)

if nargin<5, opts = struct; end
corr_type    = getf(opts,'corr_type','Spearman');
absval       = getf(opts,'absval',false);
figpos_bars  = getf(opts,'figpos_bars',[100 100 1100 420]);
labelON      = getf(opts,'labelON',true);
make_surfs   = getf(opts,'make_surfaces',false);
surf_bands   = getf(opts,'surf_bands',[]);
surf_level   = getf(opts,'surf_level','L2_both');
surf_title   = getf(opts,'surf_title','Cortical Features');
surf_cmax    = getf(opts,'surf_cmax',[]);
surf_cmap    = getf(opts,'surf_cmap',[]);
match_clim   = getf(opts,'matchClim',1);

% Settings
if isempty(surf_cmap), surf_cmap=colormap; end
labels_CMs   = stats_all.meta.label_comm_models;
labels_int   = stats_all.meta.interaction_label;
fnt          = 18;
switch lower(corr_type)
        case 'pearson'
            corropt= {'Type' 'Pearson'}; corrstr='Pearson''s r';
        otherwise
            corropt= {'Type' 'Spearman'}; corrstr='Spearman''s \rho'; 
end

if isfield(stats_all.meta,'effect_mode') 
    labels_effect = stats_all.meta.effect_mode;
else
    labels_effect = '';
end



% Available & requested levels
candidates = {'L1_route','L1_diff','L2_both'};
present = candidates(~cellfun(@(f) ~isfield(stats_all,f) || isempty(stats_all.(f)) , candidates));
levels  = getf(opts,'levels',present);
levels  = levels(ismember(levels,present));   % sanitize

assert(~isempty(levels), 'No usable levels found in stats_all.');

% Bands subset
Nfc = numel(ylbl);
bands = getf(opts,'bands',1:Nfc);
bands = bands(bands>=1 & bands<=Nfc);
assert(~isempty(bands),'No usable bands selected.');

% Size sanity (use first present level)
tmpl   = stats_all.(levels{1});
M_tmpl = size(tmpl.node_deltaR2,1);
Nnode  = size(tmpl.node_deltaR2,2);
assert(numel(sa_axis)==Nnode,'sa_axis length must equal number of nodes.');

% Decide if we are in the "M=1 everywhere" case (special bar layout)
Meach = cellfun(@(f) size(stats_all.(f).node_deltaR2,1), levels);
M_eq1 = all(Meach==1);

% Correlations
sa_axis=sa_axis(:); % corr requires matching column orientation.
allsurfh={}; % Defined even when make_surfaces=false.
out = struct; out.r = struct(); out.p = struct(); out.fig_surf = gobjects(0);
for L = 1:numel(levels)
    f = levels{L};
    S = stats_all.(f);
    [M, ~, ~] = size(S.node_deltaR2);
    r = nan(M, Nfc); p = nan(M, Nfc);
    for m = 1:M
        X = squeeze(S.node_deltaR2(m,:,:));   % Nnode x Nfc
        for ii = 1:Nfc
            if ~ismember(ii, bands), continue; end
            v = X(:,ii);
            valid = isfinite(v) & isfinite(sa_axis(:));
            if nnz(valid) < 5
                r(m,ii) = NaN; p(m,ii) = NaN;
            else
                [r(m,ii), p(m,ii)] = corr(v(valid), sa_axis(valid), corropt{:});
            end
        end
    end
    if absval, r = abs(r); end
    out.r.(f) = r; out.p.(f) = p;
end

% ---------- BAR PLOTS ----------
fig = myfig([labels_CMs ' (' labels_int ',  ' labels_effect '): \DeltaR^2 vs S-A axis'],figpos_bars); hold on
% colors: match viz_BICS_dual M=1 ribbons (lines(nlev) in level order)
if M_eq1
    nlev = numel(levels);
    cols = getf(opts,'m_eq1_colors',[]);
    if isempty(cols), cols = lines(nlev); end
    % Assemble Nfc x nlev matrix of r (M=1 squeeze)
    Rmat = nan(Nfc, nlev);
    for L = 1:nlev
        rL = out.r.(levels{L});  % 1 x Nfc
        Rmat(:,L) = rL(1,:).';
    end
    % Keep selected bands only
    Rmat = Rmat(bands,:);
    bands_lbl = ylbl(bands);
  % BAR PLOT
    bh = bar(Rmat, 'grouped'); grid on; set(gca,'XGrid','off'); box on;
    for k=1:numel(bh), set(bh(k),'FaceColor',cols(k,:)); end
    xticks(1:numel(bands)); font(fnt+2); ylim([-.99 .99]);
    for ii=1:numel(bands)-1; line([ii+.5 ii+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
    if labelON 
        xticklabels(bands_lbl); xtickangle(0);
        ylabel(sprintf('%s',corrstr));
        title('Nodal \DeltaR^2 vs S-A axis');
        legend(nice_names(levels), 'Location','northeastoutside'); 
    end
else
    % General case (M may differ by level). Plot per-level stacked subplots.
    tl = tiledlayout(numel(levels),1,'TileSpacing','compact','Padding','compact');
    cols = lines(max(Meach));
    for L = 1:numel(levels)
        nexttile(tl,L); hold on
        f  = levels{L};
        rL = out.r.(f);           % M x Nfc
        M  = size(rL,1);
        Rmat = rL(:,bands).';     % (bands) x M
        bh = bar(Rmat,'grouped'); grid on; set(gca,'XGrid','off'); box on;
        for k=1:min(M,numel(bh)), set(bh(k),'FaceColor',cols(k,:)); end
        xticks(1:numel(bands)); xticklabels(ylbl(bands)); xtickangle(0);
        ylabel(sprintf('%s',corrstr)); title(nice_names({f})); 
        % for ii=1:numel(bands)-1; line([ii+.5 ii+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end    
        if L==1
            if isfield(stats_all.(f),'model_labels')
                legend(stats_all.(f).model_labels, 'Location','northeastoutside');
            end
        end
    end
    title(tl, 'Nodal \DeltaR^2 vs S-A axis');
end
out.fig_bars = fig;

% ---------- OPTIONAL SURFACES ----------
if make_surfs
    % Choose which bands to surface
    if isempty(surf_bands), surf_bands = bands; end
    % Which level to show when M=1 case
    show_level = surf_level;
    if ~ismember(show_level, levels), show_level = levels{min( numel(levels), 1)}; end
    S = stats_all.(show_level);
    % Build ΔR² node maps matrix D (Nnode x |surf_bands|) for the first model (or average over M)
    if size(S.node_deltaR2,1)==1
        D = squeeze(S.node_deltaR2(1,:,surf_bands));      % Nnode x nb
    else
        D = squeeze(mean(S.node_deltaR2(:, :, surf_bands), 1)); % average over models
    end
    if isvector(D), D = D(:); end
    nb = size(D,2);

    if match_clim
        % Determine CMAX if not given
        if isempty(surf_cmax)
            vals = abs(D(:)); vals = vals(isfinite(vals));
            if isempty(vals), CMAX = 1; else, CMAX = max(.1,ceil(prctile(vals,95)*10)/10); end
        else
            CMAX = surf_cmax;
        end
    end

    % Plot each requested band
    out.fig_surf = gobjects(0,1); % Figure count comes from the surface helper.
    ttl = sprintf('%s (%s, %s): %s — %s', labels_CMs, labels_int, labels_effect, surf_title, pretty_level(show_level));
    Hcx = plot_conn_surf(D, pinfo, 'cortex', ylbl(surf_bands), ttl);
    if match_clim
        setsurf(Hcx, [0 CMAX], surf_cmap);
    else
        setsurf(Hcx, [], surf_cmap);
    end
    for ii = 1 : length(Hcx), out.fig_surf(ii) = Hcx{ii}.figure; end
    allsurfh=Hcx;
end

% Meta
out.meta = struct('levels',{levels}, 'bands',bands, 'corr_type',corr_type, 'absval',absval);

end

% ---- helpers ----
function v = getf(s, name, default)
if ~isstruct(s) || ~isfield(s,name) || isempty(s.(name)), v = default; else, v = s.(name); end
end

function c = nice_names(levels)
map = containers.Map( ...
    {'L1_route','L1_diff','L2_both','L3_int'}, ...
    {'L1-Route','L1-Diff','L2-Both','L3-Both+Int'} );
c = cellfun(@(f) map(f), levels, 'uni', 0);
end

function s = pretty_level(f)
switch f
    case 'L1_route', s = 'L1-Route';
    case 'L1_diff',  s = 'L1-Diff';
    case 'L2_both',  s = 'L2-Both';
    case 'L3_int',   s = 'L3-Both+Int';
    otherwise,       s = f;
end
end
