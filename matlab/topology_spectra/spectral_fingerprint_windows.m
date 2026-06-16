function out = spectral_fingerprint_windows(ALLCM, CMLBL, Vbasis, wins, varargin)
% SPECTRAL_FINGERPRINT_WINDOWS
% Spectral window fingerprinting for communication matrices (routing & diffusion).
% Supports AUTO windows via Laplacian smoothness quantiles.
%
% Inputs
%   ALLCM   : 1xJ cell, each {j} is 1xM cell of NxN matrices (models)
%   CMLBL   : 1xM cellstr of model labels (e.g., {'SPE','NE','SIE','PT','CMY','DE'})
%   Vbasis  : struct with spectral bases:
%               Vbasis.adj.V = NxKadj (Adj basis)
%               Vbasis.lap.V = NxKlap (Lap basis)
%   wins    : 1xW cell of index vectors, OR []/'auto' to compute automatically
%
% Name-Value options
%   'RefIdx'            : reference dataset (default 1)
%   'OpForModel'        : 1xM {'adj'|'lap'} per model (default {'adj',...,'lap'})
%   'DatasetLabels'     : 1xJ cellstr
%   'Cmap'              : Jx3 RGB
%   'Plot'              : 'bar'|'none' (default 'bar')
%   'RenormCovered'     : true/false (default true)
%   'DeoverlapWindows'  : true/false (default false)
%   'DE_UseDhalf'       : true/false (default true)
%   'Dhalf'             : D^{1/2} matrix for DE congruence
%   'LabelMode'         : 'on'|'off' (default 'on')
%   'FontSize'          : default 20
%   'FontName'          : default 'Calibri'
%   'AutoWindows'       : true/false (default auto if wins==[] or 'auto')
%   'Lsym'              : normalized Laplacian (required if AutoWindows)
%   'QuantileEdges'     : vector (default [0 .05 .25 .55 1])  % 4 windows
%   'StartAt'           : integer (default 2)  % exclude k=1 from windowing
%   'SplitDC'           : make a separate DC (k=1) panel
%   'ShowBothDC'        : show DC for both bases (Adj & Lap)
%
% Output: out.phi, out.phi_raw, out.coverage, out.dphi, out.wins, etc.
%
%
% 2025 Mark C Nelson (MNI)
%--------------------------------------------------------------------------


    % ---------- parse options ----------
    p = inputParser;
    addParameter(p,'SplitDC',true,@islogical);
    addParameter(p,'ShowBothDC',true,@islogical);
    addParameter(p,'RefIdx',1,@(x)isnumeric(x)&&isscalar(x));
    addParameter(p,'OpForModel',{},@(x)iscell(x));
    addParameter(p,'DatasetLabels',{},@(x)iscell(x));
    addParameter(p,'Cmap',[],@(x)isnumeric(x)&&size(x,2)==3);
    addParameter(p,'Plot','bar',@(s)ischar(s)||isstring(s));
    addParameter(p,'RenormCovered',true,@islogical);
    addParameter(p,'DeoverlapWindows',false,@islogical);
    addParameter(p,'DE_UseDhalf',true,@islogical);
    addParameter(p,'Dhalf',[],@(x)ismatrix(x) && size(x,1)==size(x,2));
    addParameter(p,'LabelMode','on',@(s)ischar(s)||isstring(s));
    addParameter(p,'FontSize',20,@(x)isnumeric(x)&&isscalar(x));
    addParameter(p,'FontName','Calibri',@(s)ischar(s)||isstring(s));
    addParameter(p,'AutoWindows',[],@(x)islogical(x) || isempty(x));
    addParameter(p,'Lsym',[],@(x)ismatrix(x) && size(x,1)==size(x,2));
    addParameter(p,'QuantileEdges',[0 .05 .25 .55 1],@(x)isnumeric(x)&&isvector(x)&&x(1)==0&&x(end)==1);
    addParameter(p,'StartAt',2,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
    parse(p,varargin{:});

    SplitDC           = p.Results.SplitDC;
    ShowBothDC        = p.Results.ShowBothDC;
    RefIdx            = p.Results.RefIdx;
    OpForModel        = p.Results.OpForModel;
    dsLabels          = p.Results.DatasetLabels;
    cmap              = p.Results.Cmap;
    plotMode          = lower(string(p.Results.Plot));
    renormCovered     = p.Results.RenormCovered;
    deoverlapWins     = p.Results.DeoverlapWindows;
    DE_UseDhalf       = p.Results.DE_UseDhalf;
    Dhalf             = p.Results.Dhalf;
    labelMode         = lower(string(p.Results.LabelMode));
    fntsz             = p.Results.FontSize;
    fntname           = char(p.Results.FontName);
    AutoWindowsOpt    = p.Results.AutoWindows;
    Lsym              = p.Results.Lsym;
    qEdges            = p.Results.QuantileEdges(:).';
    startAt           = p.Results.StartAt;

    % ---------- sizes & labels ----------
    J = numel(ALLCM);
    M = numel(ALLCM{1});
    assert(iscell(CMLBL) && numel(CMLBL)==M, 'CMLBL must be 1xM cellstr.');

    if isempty(dsLabels), dsLabels = arrayfun(@(j_) sprintf('D%d', j_), 1:J, 'UniformOutput', false); end
    if isempty(cmap), cmap = lines(J); end

    if isempty(OpForModel)
        OpForModel = repmat({'adj'}, 1, M);
        OpForModel{end} = 'lap';
    else
        assert(numel(OpForModel)==M, 'OpForModel must be 1xM.');
    end

    % ---------- bases ----------
    assert(isfield(Vbasis,'adj') && isfield(Vbasis.adj,'V'), 'Vbasis.adj.V missing');
    assert(isfield(Vbasis,'lap') && isfield(Vbasis.lap,'V'), 'Vbasis.lap.V missing');
    Vadj = Vbasis.adj.V;
    Vlap = Vbasis.lap.V;
    Kadj = size(Vadj,2);
    Klap = size(Vlap,2);

    % --- DC projectors (k=1) ---
    % Adj DC: Perron (strength/hub) mode
    v1_adj   = Vadj(:,1);
    P_dc_adj = v1_adj * v1_adj.';
    
    % Lap DC: normalized-graph constant: u ∝ D^{1/2} * 1
    if ~isempty(Dhalf)
        u = Dhalf * ones(size(Vlap,1),1);
        u = u / max(norm(u), eps);
    else
        % fallback if congruence disabled everywhere
        u = Vlap(:,1);
        u = u / max(norm(u), eps);
    end
    P_dc_lap = u * u.';

    % ---------- AUTO WINDOWS (Lap smoothness quantiles) ----------
    if isempty(AutoWindowsOpt)
        AutoWindows = (isempty(wins) || (ischar(wins) && strcmpi(wins,'auto')) || (isstring(wins) && lower(wins)=="auto"));
    else
        AutoWindows = AutoWindowsOpt;
    end

    if AutoWindows
        assert(~isempty(Lsym), 'AutoWindows requires ''Lsym'' (normalized Laplacian).');
        K = min(Klap, size(Lsym,1));                  % use first Klap columns
        Vuse = Vlap(:,1:K);
        % Smoothness scores: S_k = v_k^T L v_k
        S = sum((Lsym*Vuse).*Vuse,1);                 % 1 x K
        % Sort ascending (low S = global)
        [~, ord] = sort(S, 'ascend');
        cutIdx = unique(max(1, min(K, round(qEdges*K))));
        winsLap = cell(1, numel(cutIdx)-1);
        for i=1:numel(cutIdx)-1
            slab = ord(cutIdx(i)+1 : cutIdx(i+1));    % indices in 1..K
            % Enforce StartAt (exclude k < startAt from windowing)
            slab = slab(slab >= startAt);
            winsLap{i} = sort(slab);
        end
        % Clean empty windows
        winsLap = winsLap(~cellfun(@isempty,winsLap));
        % Apply these SAME index sets to adjacency as well (scale-consistent bins)
        winsAdj = winsLap;
        wins = winsLap;  % stored for reporting
    else
        % user-supplied windows (apply same sets to both bases)
        assert(iscell(wins) && ~isempty(wins), 'wins must be a non-empty cell array, or use AutoWindows.');
        winsAdj = wins; winsLap = wins;
        % Enforce StartAt for both
        for w=1:numel(winsAdj), winsAdj{w} = winsAdj{w}(winsAdj{w}>=startAt); end
        for w=1:numel(winsLap), winsLap{w} = winsLap{w}(winsLap{w}>=startAt); end
    end
    W = numel(winsLap);

    % ---------- optional de-overlap for fractions ----------
    wins_adj_uniq = ifelse_deoverlap(winsAdj, Kadj, deoverlapWins);
    wins_lap_uniq = ifelse_deoverlap(winsLap, Klap, deoverlapWins);

    % ---------- projectors ----------
    P_adj      = cell(1,W); P_lap      = cell(1,W);
    P_adj_uniq = cell(1,W); P_lap_uniq = cell(1,W);
    for w=1:W
        ia = clamp_idx(winsAdj{w}, Kadj);
        il = clamp_idx(winsLap{w}, Klap);
        P_adj{w} = Vadj(:,ia)*Vadj(:,ia).';
        P_lap{w} = Vlap(:,il)*Vlap(:,il).';

        iau = clamp_idx(wins_adj_uniq{w}, Kadj);
        ilu = clamp_idx(wins_lap_uniq{w}, Klap);
        P_adj_uniq{w} = Vadj(:,iau)*Vadj(:,iau).';
        P_lap_uniq{w} = Vlap(:,ilu)*Vlap(:,ilu).';
    end
  % UNION (no double counting) — force k=1 in union
    idxUnion_adj = clamp_idx(unique([1, winsAdj{:}]), Kadj);
    idxUnion_lap = clamp_idx(unique([1, winsLap{:}]), Klap);  
    P_adj_union = Vadj(:,idxUnion_adj)*Vadj(:,idxUnion_adj).';
    P_lap_union = Vlap(:,idxUnion_lap)*Vlap(:,idxUnion_lap).';

    % ---------- storage ----------
    phi_raw  = nan(W,J,M);
    phi      = nan(W,J,M);
    coverage = nan(1,J,M);
    dphi     = nan(W,J,M);

    % ---------- core ----------
    for j=1:J
        for m=1:M
            Mij = ALLCM{j}{m};
            Mij = symm_nan0(Mij);

            op = lower(OpForModel{m});
            if strcmp(op,'lap') && DE_UseDhalf
                if isempty(Dhalf), error('Provide Dhalf (D^{1/2}) when DE_UseDhalf=true.'); end
                Mij = Dhalf * Mij * Dhalf;     % map to Lsym geometry
            end
            totalE = fro2_nan(Mij); if totalE<=0 || ~isfinite(totalE), continue; end

            % fractions (optionally de-overlapped)
            E_sum = 0;
            for w=1:W
                Pw = iff(strcmp(op,'adj'), P_adj_uniq{w}, P_lap_uniq{w});
                eW = fro2_nan(Pw*Mij*Pw);
                phi_raw(w,j,m) = eW / totalE;
                E_sum = E_sum + eW;
            end

            % coverage via UNION projector
            Punion = iff(strcmp(op,'adj'), P_adj_union, P_lap_union);
            coverage(1,j,m) = fro2_nan(Punion*Mij*Punion) / totalE;

            % renormalize within covered spectrum
            if renormCovered && E_sum>0
                phi(:,j,m) = phi_raw(:,j,m) * (totalE / E_sum);
            else
                phi(:,j,m) = phi_raw(:,j,m);
            end
        end
    end

  % ------- Compute DC fractions per dataset & model ----------
    DC_adj = nan(J,M);
    DC_lap = nan(J,M);
    
    for j = 1:J
      for m = 1:M
        Mij = symm_nan0(ALLCM{j}{m});
    
        % Adj DC fraction (Euclidean)
        EtotA = fro2_nan(Mij);
        if EtotA > 0
            DC_adj(j,m) = fro2_nan(P_dc_adj * Mij * P_dc_adj) / EtotA;
        end
    
        % Lap DC fraction (always in Lsym geometry)
        Mlap  = Mij;
        if DE_UseDhalf && ~isempty(Dhalf)
            Mlap = Dhalf * Mij * Dhalf;
        end
        EtotL = fro2_nan(Mlap);
        if EtotL > 0
            DC_lap(j,m) = fro2_nan(P_dc_lap * Mlap * P_dc_lap) / EtotL;
        end
      end
    end

    % ---------- deltas vs reference ----------
    for m=1:M
        refcol = phi(:,RefIdx,m);
        for j=1:J, dphi(:,j,m) = phi(:,j,m) - refcol; end
    end

    % ---------- pack ----------
    out = struct();
    out.phi          = phi;
    out.phi_raw      = phi_raw;
    out.coverage     = coverage;
    out.dphi         = dphi;
    out.wins         = wins;               % reported (pre de-overlap)
    out.OpForModel   = OpForModel;
    out.DatasetLabels= dsLabels;
    out.ModelLabels  = CMLBL;
    out.fig          = [];

    % ---------- plot ----------
    % barwidth=0.86;
    barwidth=0.75;
    if plotMode=="bar"
        try
            out.fig = myfig('Spectral fingerprints — window fractions & deltas', [80 80 1800 820]);
        catch
            out.fig = figure('Position',[80 80 1800 820],'Color','w');
        end
        T = tiledlayout(2, numel(wins), 'Padding','compact','TileSpacing','compact');

        % handles
        axTop = gobjects(1,numel(wins));
        axBot = gobjects(1,numel(wins));

        % --- Top row: FRACTIONS
        for w=1:numel(wins)
            axTop(w) = nexttile; hold on;
            Mmat  = squeeze(phi(w,:,:)).'; Mplot = Mmat; Mplot(~isfinite(Mplot))=0;
            B = bar(Mplot,'grouped','BarWidth',barwidth);
            for jdx=1:J, B(jdx).FaceColor='flat'; B(jdx).CData = repmat(cmap(jdx,:), size(Mplot,1),1); end
            if renormCovered, ylim([0 1]); yline(1,'k:'); else, ylim([0 1]); end
            xlim([0.5,size(Mplot,1)+0.5]); yl=ylim; for s=1:(size(Mplot,1)-1), line([s+0.5 s+0.5],yl,'Color','k','LineStyle','-','LineWidth',0.5,'HitTest','off'); end
            if strcmpi(labelMode,'on')
                if w==1, ylabel('Fraction in window'); end
                set(axTop(w),'XTick',1:size(Mplot,1),'XTickLabel',[]);
                title(axTop(w), win_label(wins{w}));
                if w==1, legend(B, dsLabels,'Location','northoutside','Orientation','horizontal'); end
            else
                ylabel(axTop(w),''); set(axTop(w),'XTickLabel',[]); title(axTop(w),''); legend off;
                if w~=1, set(axTop(w),'YTickLabel',''); end
            end
            set(axTop(w),'XGrid','off','YGrid','on'); box on; font(fntsz,fntname); hold off;
        end

        % --- Bottom row: DELTAS
        ymax=0; for w=1:numel(wins), X=squeeze(dphi(w,:,:)); ymax=max(ymax, max(abs(X(:)),[],'omitnan')); end
        ylims = (ymax==0)*[-0.1 0.1] + (ymax>0)*[-1.05*ymax 1.05*ymax];
        for w=1:numel(wins)
            axBot(w) = nexttile; hold on;
            Mmat  = squeeze(dphi(w,:,:)).'; Mplot = Mmat; Mplot(~isfinite(Mplot))=0;
            B = bar(Mplot,'grouped','BarWidth',barwidth);
            for jdx=1:J, B(jdx).FaceColor='flat'; B(jdx).CData = repmat(cmap(jdx,:), size(Mplot,1),1); end
            yline(0,'k-'); ylim(ylims);
            xlim([0.5,size(Mplot,1)+0.5]); yl=ylim; for s=1:(size(Mplot,1)-1), line([s+0.5 s+0.5],yl,'Color','k','LineStyle','-','LineWidth',0.5,'HitTest','off'); end
            if strcmpi(labelMode,'on')
                if w==1, ylabel('\Delta fraction (vs ref)'); end
                set(axBot(w),'XTick',1:size(Mplot,1),'XTickLabel',shorten_labels(CMLBL));
            else
                ylabel(axBot(w),''); set(axBot(w),'XTickLabel',[]);
                if w~=1, set(axBot(w),'YTickLabel',''); end
            end
            set(axBot(w),'XGrid','off','YGrid','on'); box on; font(fntsz,fntname); hold off;
        end

        if strcmpi(labelMode,'on')
            covTxt = coverage_summary(coverage, dsLabels, CMLBL);
            title(T, sprintf('Spectral fingerprints (autoWins=%d, renorm=%d, deoverlap=%d)\n%s', ...
                AutoWindows, renormCovered, deoverlapWins, covTxt));
            font(fntsz,fntname);
        end
    end

  % ------------ plot DC (k=1) -----------------
    if SplitDC
        try
            figDC = myfig('DC (k=1) contributions — Adj & Lap geometries', [120 80 1400 520]);
        catch
            figDC = figure('Position',[120 80 1400 520],'Color','w');
        end
        Tdc = tiledlayout(1, ShowBothDC + 1, 'Padding','compact','TileSpacing','compact');
    
        % Panel 1: Adj DC fractions (per model, grouped by dataset colors)
        nexttile; hold on;
        Mplot = DC_adj.'; Mplot(~isfinite(Mplot)) = 0;   % M x J
        B = bar(Mplot, 'grouped', 'BarWidth', barwidth);
        for jdx=1:J, B(jdx).FaceColor='flat'; B(jdx).CData = repmat(cmap(jdx,:), size(Mplot,1), 1); end
        ylim([0 1]); xlim([0.5 size(Mplot,1)+0.5]); yl=ylim;
        for s=1:(size(Mplot,1)-1), line([s+0.5 s+0.5], yl, 'Color','k','LineWidth',0.5); end
        if strcmpi(labelMode,'on')
            ylabel('Fraction in k=1 (Adj geometry)'); set(gca,'XTick',1:size(Mplot,1),'XTickLabel',shorten_labels(CMLBL)); 
            legend(B, dsLabels, 'Location','northoutside','Orientation','horizontal');
            title('Adjacency DC (v_1) contribution');
        else
            set(gca,'XTickLabel',''); title(''); legend off;
        end
        set(gca,'XGrid','off','YGrid','on'); box on; font(fntsz,fntname); hold off;
    
        if ShowBothDC
            % Panel 2: Lap DC fractions (per model, grouped by dataset colors)
            nexttile; hold on;
            Mplot = DC_lap.'; Mplot(~isfinite(Mplot)) = 0;   % M x J
            B = bar(Mplot, 'grouped', 'BarWidth', barwidth);
            for jdx=1:J, B(jdx).FaceColor='flat'; B(jdx).CData = repmat(cmap(jdx,:), size(Mplot,1), 1); end
            ylim([0 1]); xlim([0.5 size(Mplot,1)+0.5]); yl=ylim;
            for s=1:(size(Mplot,1)-1), line([s+0.5 s+0.5], yl, 'Color','k','LineWidth',0.5); end
            if strcmpi(labelMode,'on')
                ylabel('Fraction in k=1 (Lap geometry)'); set(gca,'XTick',1:size(Mplot,1),'XTickLabel',shorten_labels(CMLBL));
                title('Laplacian DC (constant) contribution');
            else
                set(gca,'XTickLabel',''); title(''); set(gca,'YTickLabel','');
            end
            set(gca,'XGrid','off','YGrid','on'); box on; font(fntsz,fntname); hold off;
        end
    end

% -------------------------------------------------------------------------
end

% ----------------- helpers -----------------
function sets_uniq = ifelse_deoverlap(sets, K, doit)
    if ~doit, sets_uniq = sets; return; end
    seen = false(1,K); sets_uniq = sets;
    for w=1:numel(sets)
        idx = clamp_idx(sets{w},K); keep = ~seen(idx);
        sets_uniq{w} = idx(keep); seen(idx)=true;
    end
end

function idx = clamp_idx(idx, K)
    idx = unique(idx(:))'; idx = idx(idx>=1 & idx<=K);
end

function X = symm_nan0(X)
    X = (X+X')/2; X(~isfinite(X)) = 0;
end

function E = fro2_nan(M)
    X = M(:); X(~isfinite(X)) = 0; E = sum(X.^2);
end

function L = shorten_labels(Lin)
    L = Lin;
    for i=1:numel(L)
        s = string(L{i});
        if strlength(s) > 16, L{i} = char(extractBefore(s,16)+"…"); end
    end
end

function txt = coverage_summary(cov, dsLabels, mdlLabels)
    C = squeeze(cov); med_by_m = median(C,1,'omitnan');
    parts = strings(1,numel(mdlLabels));
    for m=1:numel(mdlLabels), parts(m) = sprintf('%s medCov=%.2f', mdlLabels{m}, med_by_m(m)); end
    txt = strjoin(parts,' | ');
end

function s = win_label(idx)
    if isempty(idx), s = 'empty'; else, s = sprintf('k %d–%d', min(idx), max(idx)); end
end

function y = iff(cond, a, b)
    if cond, y=a; else, y=b; end
end
