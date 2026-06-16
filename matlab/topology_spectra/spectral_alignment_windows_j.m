function out = spectral_alignment_windows_j(A_ref,A_x_stack,labels,str,opts)
% SPECTRAL_ALIGNMENT_WINDOWS_J
% Compare a reference network against J target networks using windowed
% subspace alignment for adjacency & Laplacian operators.
%
% Global banded matching across k, with windowed summaries.
%
% Inputs
%   A_ref      : NxN reference adjacency (weighted, symmetric recommended)
%   A_x_stack  : NxNxJ stack of target adjacencies (each compared to A_ref)
%   labels     : 1xJ cell of data labels
%   str        : char label for figure titles
%   opts       : (optional) struct with fields:
%                  .kmax        (default 50)
%                  .normalize   (default true)  % Ahat = D^-1/2 A D^-1/2; Lsym = I - Ahat
%                  .wins        (default {2:5, 5:10, 10:15, 15:25, 25:50})
%                  .doSignAlign (default true)
%                  .coords      (Nx3) node coordinates for optional maps
%                  .SA_axis     (Nx1) for SA correlation metric
%                  .modules     (Nx1 int) module labels for module R^2
%                  .Wspatial    (NxN) spatial weight matrix for Moran's I
%                  .plotMaps    (default false)
%                  .cmap        (Jx3) colors for datasets (grouped bars)
%                  .matchMode   'global-banded' (default), or 'window-local'
%                  .band        integer band half-width around diagonal for global matching (default: max(3, round(0.06*kmax)))
%                  .LabelMode   'on'|'off' (default 'on')  % hide/show titles, labels, ticklabels, legends
%                  .KPerPage    positive int (default 100) % split per-k plots into pages if K > KPerPage
%
% Displays
%   Fig 1: Windowed mean cosine(bar) for Adj & Lap (two subplots)
%          X groups = windows; within each group J bars colored by cmap.
%   Fig 2+: Per-eigenvector grouped bars of matched cosine (Adj & Lap)
%          Automatically paged if number of k exceeds KPerPage.
%
% Output
%   out.adj / out.lap:
%       .Vref, .Vx{j}, .eval_ref, .eval_x{j}
%       .subspace_aff{win}(j) : mean cosine of principal angles (window, dataset)
%       .match{win}(j).idxRef, .idxX, .cosMatched (per-k matched cosine)
%       .topology{win}(j,k).metrics : optional topology metrics
%   out.figs.subspace, out.figs.per_k (vector of handles if paged)
%   out.labels : dataset descriptions
%   out.cmap   : colormap for target datasets
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
% -------------------------------------------------------------------------

    if nargin < 5, opts = struct(); end
    kmax        = getOpt(opts,'kmax',50);
    normalize   = getOpt(opts,'normalize',true);
    wins        = getOpt(opts,'wins',{2:5, 5:10, 10:15, 15:25, 25:50});
    doSignAlign = getOpt(opts,'doSignAlign',true);
    coords      = getOpt(opts,'coords',[]);
    SA_axis     = getOpt(opts,'SA_axis',[]);
    modules     = getOpt(opts,'modules',[]);
    Wspatial    = getOpt(opts,'Wspatial',[]);
    plotMaps    = getOpt(opts,'plotMaps',false);
    cmap        = getOpt(opts,'cmap',[]);
    matchMode   = getOpt(opts,'matchMode','global-banded');
    bandOpt     = getOpt(opts,'band',[]);
    labelMode   = getOpt(opts,'LabelMode','on');
    showLabels  = strcmpi(labelMode,'on');
    KPerPage    = getOpt(opts,'KPerPage',100);

    [N1,N2,J] = size(A_x_stack);
    assert(N1==N2 && all(size(A_ref)==[N1 N2]), 'Dimensions must be NxN and consistent.');
    if ~isempty(cmap), assert(size(cmap,1)==J && size(cmap,2)==3, 'cmap must be Jx3.'); end

    % Settings
      fntsz=20; fntname='Calibri';

    % ---- Build operators
    if normalize
        [Ahat_ref, Lsym_ref] = normalize_operators(A_ref);
        Ahat_x  = cell(1,J); Lsym_x = cell(1,J);
        for j=1:J, [Ahat_x{j}, Lsym_x{j}] = normalize_operators(A_x_stack(:,:,j)); end
        adj_label = 'Adj (norm)'; lap_label = 'Lap (norm)';
    else
        Ahat_ref = symm(A_ref);
        Lsym_ref = diag(sum(A_ref,2)) - Ahat_ref;
        Ahat_x  = cell(1,J); Lsym_x = cell(1,J);
        for j=1:J
            Ax = symm(A_x_stack(:,:,j));
            Ahat_x{j} = Ax;
            Lsym_x{j} = diag(sum(Ax,2)) - Ax;
        end
        adj_label = 'Adj'; lap_label = 'Lap';
    end

    % ---- Eigendecompositions
    [VA_ref, DA_ref] = eigs(Ahat_ref, kmax, 'largestabs');  evalA_ref = diag(DA_ref);
    [VL_ref, DL_ref] = safe_eigs_laplacian(Lsym_ref, kmax); evalL_ref = diag(DL_ref);

    VA_x = cell(1,J); evalA_x = cell(1,J);
    VL_x = cell(1,J); evalL_x = cell(1,J);
    for j=1:J
        [VAj, DAj] = eigs(Ahat_x{j}, kmax, 'largestabs');  VA_x{j} = VAj; evalA_x{j} = diag(DAj);
        [VLj, DLj] = safe_eigs_laplacian(Lsym_x{j}, kmax); VL_x{j} = VLj; evalL_x{j} = diag(DLj);
    end

    % ---- Compute windowed summaries using GLOBAL BANDED matching
    bandA = bandOpt; if isempty(bandA), bandA = max(3, round(0.06*kmax)); end
    bandL = bandA;

    adj = compute_windows_global(VA_ref, VA_x, Lsym_ref, wins, doSignAlign, SA_axis, modules, Wspatial, coords, plotMaps, sprintf('%s | %s',str,adj_label), matchMode, bandA);
    lap = compute_windows_global(VL_ref, VL_x, Lsym_ref, wins, doSignAlign, SA_axis, modules, Wspatial, coords, plotMaps, sprintf('%s | %s',str,lap_label), matchMode, bandL);

    % ================== FIGURE 1: Windowed mean cosine ====================
    fig1 = myfig(['Spectral windows: ' str ' (mean cos (principal angles))'],[300 250 1020 460]);
    t = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; 
        plot_window_bars(adj, wins, cmap, labels, showLabels);
        title_if(showLabels, adj_label);
        grid on; set(gca,'XGrid','off'); box on; font(fntsz,fntname);
    nexttile; 
        plot_window_bars(lap, wins, cmap, labels, showLabels);
        title_if(showLabels, lap_label);
        grid on; set(gca,'XGrid','off'); box on; font(fntsz,fntname);

    if showLabels
        ylabel(t,'Mean cosine','FontName',fntname,'FontSize',fntsz);
    end

    % ========== FIGURE 2+: Matched alignment per-k (with paging) ==========
    K = size(adj.cos_k,1);                       % number of eigenvectors
    nPages = max(1, ceil(K / KPerPage));
    per_k_figs = gobjects(nPages,1);

    for pg = 1:nPages
        k1 = (pg-1)*KPerPage + 1;
        k2 = min(pg*KPerPage, K);
        kShow = k1:k2;

        per_k_figs(pg) = myfig(sprintf('Matched alignment per-k: %s  (k %d–%d) [page %d/%d]', str, k1, k2, pg, nPages), [300 740 1020 520]);
        t2 = tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

        nexttile; 
            plot_perk_grouped(adj, cmap, kShow, showLabels);
            title_if(showLabels, [adj_label ' — matched |cos| per k']);
            grid on; set(gca,'XGrid','off'); box on; ylim([0 1]); font(fntsz,fntname);

        nexttile; 
            plot_perk_grouped(lap, cmap, kShow, showLabels);
            title_if(showLabels, [lap_label ' — matched |cos| per k']);
            grid on; set(gca,'XGrid','off'); box on; ylim([0 1]); font(fntsz,fntname);

        if showLabels
            xlabel(t2,'Eigenvector index (k)'); 
            ylabel(t2,'|cos| matched','FontName',fntname,'FontSize',fntsz);
        end
    end

    % ---- Output
    out = struct();
    out.adj = struct('Vref',VA_ref,'Vx',{VA_x},'eval_ref',evalA_ref,'eval_x',{evalA_x}, ...
                     'subspace_aff',{adj.subspace_aff}, 'match',{adj.match}, ...
                     'topology',{adj.topology}, 'perm',{adj.perm});
    out.lap = struct('Vref',VL_ref,'Vx',{VL_x},'eval_ref',evalL_ref,'eval_x',{evalL_x}, ...
                     'subspace_aff',{lap.subspace_aff}, 'match',{lap.match}, ...
                     'topology',{lap.topology}, 'perm',{lap.perm});
    out.figs = struct('subspace',fig1,'per_k',per_k_figs);
    out.labels=labels; out.cmap=cmap;
end

% -------------------------------------------------------------------------

% ===================== GLOBAL BANDED MATCHING WORKFLOW =====================
function S = compute_windows_global(Vref, Vx_cell, Lref, wins, doSignAlign, SA_axis, modules, Wspatial, coords, doMaps, ttl, matchMode, band)
    J = numel(Vx_cell);
    K = size(Vref,2);
    nW = numel(wins);

    % --- Sign-align to reference (optional)
    if doSignAlign
        for j=1:J, Vx_cell{j} = sign_align(Vref, Vx_cell{j}); end
    end

    % --- Build global permutations
    perm = cell(1,J);
    for j=1:J
        switch lower(matchMode)
            case 'global-banded'
                perm{j} = global_banded_match(Vref, Vx_cell{j}, band);
            case 'window-local'
                perm{j} = (1:K).';
            otherwise
                error('Unknown matchMode: %s', matchMode);
        end
    end

    % --- Permute each target basis globally to ref order 
    Vx_perm = cell(1,J);
    for j=1:J
        pj = perm{j}(:);                 % ensure column vector, length K
        if any(pj==0)
            z = find(pj==0);             % positions that were unmatched
            pj(z) = z;                   % map those ref indices to themselves
        end
        Vx_perm{j} = Vx_cell{j}(:, pj);
    end

    % --- Baseline diagonal cos (before/after perm doesn’t matter much here)
    cos_k = zeros(K,J);
    for j=1:J
        C = abs(Vref' * Vx_perm{j});
        cos_k(:,j) = diag(C);
    end

    % --- Windowed summaries (subspace + per-k matched using global perm) 
    subspace_aff = cell(1,nW);
    match = cell(1,nW);
    topology = cell(1,nW);

    for w=1:nW
        idx = wins{w};
        % Subspace affinity (use permuted Vx to align identities across windows)
        aff = zeros(1,J);
        for j=1:J
            aff(j) = subspace_affinity(Vref(:,idx), Vx_perm{j}(:,idx));
        end
        subspace_aff{w} = aff;

        % Matched pairs: each ref k maps to target perm(k)
        mcell = cell(1,J);
        topocell = cell(1,J);
        for j=1:J
            idxRef = idx(:)';                          % window ref ks
            idxX   = perm{j}(idxRef);                 % their global matches
            idxX(idxX==0) = idxRef(idxX==0);          % defensive
            cosMatched = zeros(1,numel(idxRef));
            T = repmat(struct('metrics',[]), 1, numel(idxRef));
            for t=1:numel(idxRef)
                r = idxRef(t); x = idxX(t);
                vr = Vref(:,r); vx = Vx_cell{j}(:,x); % use original Vx for metrics
                cosMatched(t) = abs(vr' * vx);
                T(t).metrics = characterize_topology(vr, vx, Lref, SA_axis, modules, Wspatial);
                % optional maps omitted here (they reference outer-scope font vars)
            end
            % local positions
            locMap = containers.Map(num2cell(idx), num2cell(1:numel(idx)));
            idxRef_local = arrayfun(@(k) locMap(k), idxRef);
            idxX_local = nan(size(idxX));
            inWin = ismember(idxX, idx);
            idxX_local(inWin) = arrayfun(@(k) locMap(k), idxX(inWin));

            mcell{j} = struct('idxRef',idxRef, 'idxX',idxX, ...
                              'idxRef_local', idxRef_local, 'idxX_local', idxX_local, ...
                              'cosMatched',cosMatched);
            topocell{j} = T;
        end
        match{w} = mcell;
        topology{w} = topocell;
    end

    S = struct('subspace_aff',{subspace_aff}, 'match',{match}, 'topology',{topology}, ...
               'cos_k',cos_k, 'perm',{perm});
end

% -------------------------------------------------------------------------

% ===================== Global banded Hungarian  ==================
function perm = global_banded_match(Vref, Vx, band)
% Returns a column vector perm of length K: ref k -> target perm(k)
    K = size(Vref,2);
    S = abs(Vref' * Vx);                 % K x K similarity
    cost = 1 - S;
    [I,J] = ndgrid(1:K,1:K);
    cost(abs(I-J) > band) = 1e6;         % big penalty outside band
    % Prefer matchpairs if available; else use munkres
    if exist('matchpairs','file') == 2
        pairs = matchpairs(cost, 1e9);   % very high unmatched cost threshold
        perm = zeros(K,1);
        perm(pairs(:,1)) = pairs(:,2);
    else
        assign = munkres(cost);          % 1xK or Kx1 vector: assign(i)=j
        perm = assign(:);
    end
    % If any zeros slipped through (shouldn’t), greedily fill with best remaining:
    if any(perm==0)
        freeR = find(perm==0);
        usedC = perm(perm>0);
        avail = setdiff(1:K, usedC);
        for r = freeR(:)'
            [~,ix] = max(S(r,avail));
            perm(r) = avail(ix);
            avail(ix) = [];
        end
    end
end

% ======================= Plot helpers  ============
function plot_window_bars(S, wins, cmap, labels, showLabels)
    nW = numel(wins);
    J = numel(S.subspace_aff{1});
    M = zeros(nW, J);
    for w=1:nW, M(w,:) = S.subspace_aff{w}; end
    B = bar(M, 'grouped', 'BarWidth', 0.75);
    if ~isempty(cmap)
        for j=1:J
            B(j).FaceColor = 'flat';
            B(j).CData = repmat(cmap(j,:), nW, 1);
        end
    end
    lbls = arrayfun(@(i) sprintf('%d-%d', wins{i}(1), wins{i}(end)), 1:numel(wins), 'UniformOutput', false);
    if showLabels
        set(gca,'XTickLabel',lbls);
        if ~isempty(labels), legend(labels,'Location','southoutside','Orientation','horizontal'); end
    else
        set(gca,'XTickLabel',[]);
        legend off;
        title('');
        xlabel(''); ylabel('');
    end
    ylim([0 1]);
end

function plot_perk_grouped(S, cmap, kShow, showLabels)
    if nargin<3 || isempty(kShow), kShow = []; end
    if nargin<4, showLabels = true; end
    K = size(S.cos_k,1); J = size(S.cos_k,2);
    M = S.cos_k; % baseline diag cos
    % overwrite with matched global cos where available (by ref k)
    nW = numel(S.match);
    for w=1:nW
        for j=1:J
            m = S.match{w}{j};
            for t=1:numel(m.idxRef)
                kr = m.idxRef(t);
                M(kr,j) = m.cosMatched(t);
            end
        end
    end
    if ~isempty(kShow)
        kShow = kShow(:)'; kShow = kShow(kShow>=1 & kShow<=K);
        Mplot = M(kShow,:); xvals = kShow;
    else
        Mplot = M; xvals = 1:K;
    end
    B = bar(Mplot, 'grouped', 'BarWidth', 0.92);
    if ~isempty(cmap)
        for j=1:J, B(j).FaceColor = 'flat'; B(j).CData = repmat(cmap(j,:), size(Mplot,1), 1); end
    end
    xlim([0 size(Mplot,1)+1]); ylim([0 1]);
    set(gca,'XTick',1:numel(xvals));
    if showLabels
        set(gca,'XTickLabel',xvals);
    else
        set(gca,'XTickLabel',[]);
        legend off; title(''); xlabel(''); ylabel('');
    end
end

% ======================= Numerics & misc helpers ===========================
function [V,D] = safe_eigs_laplacian(L, k)
    n = size(L,1);
    k = min(k, n-1);
    opts = struct('issym', true, 'isreal', true);
    try
        [V,D] = eigs(L, k, 'smallestabs', opts);
    catch
        s = 1e-6 * max(1, full(trace(L))/n);
        [V,D] = eigs(L, k, s, opts);
    end
end
function [Ahat, Lsym] = normalize_operators(A)
    A = symm(A);
    d = sum(A,2); d(d<=0) = eps;
    Dmh = spdiags(1./sqrt(d), 0, size(A,1), size(A,2));
    Ahat = Dmh * A * Dmh;
    Lsym = speye(size(A)) - Ahat;
end
function X = symm(X), X = (X+X')/2; end
function v = getOpt(s,f,d), if isfield(s,f), v=s.(f); else, v=d; end, end
function Vx = sign_align(Vref, Vx)
    C = Vref' * Vx;
    sgn = sign(diag(C)); sgn(sgn==0)=1;
    Vx = Vx .* sgn';
end
function aff = subspace_affinity(U, V)
    [~,S,~] = svd(U' * V, 'econ'); s = diag(S);
    aff = mean(min(max(s,0),1));
end
function metrics = characterize_topology(vr, vx, Lref, SA_axis, modules, Wspatial)
    metrics = struct();
    if ~isempty(SA_axis)
        metrics.SAr_ref = corr(vr, SA_axis, 'type','Spearman','rows','pairwise');
        metrics.SAr_x   = corr(vx, SA_axis, 'type','Spearman','rows','pairwise');
    end
    metrics.smooth_ref = vr' * Lref * vr;
    metrics.smooth_x   = vx' * Lref * vx;
    metrics.PR_ref = (sum(vr.^2)^2) / sum(vr.^4);
    metrics.PR_x   = (sum(vx.^2)^2) / sum(vx.^4);
    if ~isempty(modules)
        metrics.modR2_ref = module_R2(vr, modules);
        metrics.modR2_x   = module_R2(vx, modules);
    end
    if ~isempty(Wspatial)
        metrics.MoranI_ref = moransI(vr, Wspatial);
        metrics.MoranI_x   = moransI(vx, Wspatial);
    end
end
function R2 = module_R2(v, mods)
    v = v(:); u = unique(mods(:)); X = zeros(numel(v), numel(u));
    for i=1:numel(u), X(:,i) = double(mods==u(i)); end
    b = X \ v; vhat = X*b;
    ssr = sum((vhat-mean(v)).^2); sst = sum((v-mean(v)).^2);
    R2 = max(0, min(1, ssr / max(sst, eps)));
end
function I = moransI(v, W)
    v = v(:); v0 = v - mean(v);
    num = v0' * W * v0; den = v0' * v0;
    I = num / max(den, eps);
end

% -------- label helpers (local) --------
function title_if(show, txt), if show, title(txt); else, title(''); end, end
function xlabel_if(show, txt), if show, xlabel(txt); else, xlabel(''); end, end
function ylabel_if(show, txt), if show, ylabel(txt); else, ylabel(''); end, end
function legend_if(show, varargin)
    if show, legend(varargin{:}); else, legend off; end
end
function xticklabels_if(show, vals)
    if show, set(gca,'XTickLabel',vals); else, set(gca,'XTickLabel',[]); end
end
