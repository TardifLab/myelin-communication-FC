function out = spectral_alignment_windows_robust(A_ref,A_x_stack,labels,str,opts)
% SPECTRAL_ALIGNMENT_WINDOWS_J (robust + per-k deltas; with window clipping)
% Compare a reference network against J target networks using windowed
% subspace alignment for adjacency & Laplacian operators, with:
%   • Consistent Laplacian (Lsym = I - D^-1/2 A D^-1/2)
%   • Optional DC-mode drop for Lap
%   • Global-banded matching (Hungarian), plus window summaries
%   • Procrustes alignment within windows + per-mode sign disambiguation
%   • Projector-surrogate logic retained via subspace affinity; per-mode deltas stored
%   • windows are clipped to in-range [1..K] per-space (Adj/Lap)
%   • matching band capped to K-1 per-space (prevents out-of-range)
%
%
% 2025 Mark C Nelson (MNI)
% -------------------------------------------------------------------------

    if nargin < 5, opts = struct(); end
    kmax        = getOpt(opts,'kmax',50);
    normalize   = getOpt(opts,'normalize',true);
    wins_in     = getOpt(opts,'wins',{2:5, 5:10, 10:15, 15:25, 25:50});
    doSignAlign = getOpt(opts,'doSignAlign',true);
    coords      = getOpt(opts,'coords',[]);
    SA_axis     = getOpt(opts,'SA_axis',[]);
    modules     = getOpt(opts,'modules',[]);
    Wspatial    = getOpt(opts,'Wspatial',[]);
    plotMaps    = getOpt(opts,'plotMaps',false); %#ok<NASGU>
    cmap        = getOpt(opts,'cmap',[]);
    matchMode   = getOpt(opts,'matchMode','global-banded');
    bandOpt     = getOpt(opts,'band',[]);
    labelMode   = getOpt(opts,'LabelMode','on');
    showLabels  = strcmpi(labelMode,'on');
    KPerPage    = getOpt(opts,'KPerPage',100);
    dropDC      = getOpt(opts,'dropDC',true);
    SArType     = lower(getOpt(opts,'SArType','pearson'));

    usefont='Helvetica'; fntsz=18;

    [N1,N2,J] = size(A_x_stack);
    assert(N1==N2 && all(size(A_ref)==[N1 N2]), 'Dimensions must be NxN and consistent.');
    if ~isempty(cmap), assert(size(cmap,1)==J && size(cmap,2)==3, 'cmap must be Jx3.'); end
    if isempty(SA_axis), warning('opts.SA_axis is empty: SAr-related deltas will be NaN.'); end

    % ---- Build operators (Lsym consistency) --------------------------------
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

    % ---- Eigendecompositions ------------------------------------------------
    [VA_ref, DA_ref] = eigs(Ahat_ref, kmax, 'largestabs');  evalA_ref = diag(DA_ref);
    [VL_ref, DL_ref] = safe_eigs_laplacian(Lsym_ref, kmax); evalL_ref = diag(DL_ref);
    if dropDC, [VL_ref, evalL_ref] = drop_dc_mode(VL_ref, evalL_ref); end

    VA_x = cell(1,J); evalA_x = cell(1,J);
    VL_x = cell(1,J); evalL_x = cell(1,J);
    for j=1:J
        [VAj, DAj] = eigs(Ahat_x{j}, kmax, 'largestabs');  VA_x{j} = VAj; evalA_x{j} = diag(DAj);
        [VLj, DLj] = safe_eigs_laplacian(Lsym_x{j}, kmax); VL_x{j} = VLj; evalL_x{j} = diag(DLj);
        if dropDC, [VL_x{j}, evalL_x{j}] = drop_dc_mode(VL_x{j}, evalL_x{j}); end
    end

    % Align K across datasets (after DC-drop)
    kA = min([size(VA_ref,2) cellfun(@(V) size(V,2), VA_x)]);
    kL = min([size(VL_ref,2) cellfun(@(V) size(V,2), VL_x)]);
    VA_ref = VA_ref(:,1:kA);    VA_x = cellfun(@(V) V(:,1:kA), VA_x, 'uni',0);
    VL_ref = VL_ref(:,1:kL);    VL_x = cellfun(@(V) V(:,1:kL), VL_x, 'uni',0);

    % Clip windows per space to prevent out-of-range
    winsA = clip_windows(wins_in, kA);
    winsL = clip_windows(wins_in, kL);

    % Optional diagnostics when clipping occurred
    if any(~cellfun(@(a,b) isequal(a,b), wins_in, winsA))
        warning('[ADJ] Some windows were clipped to K=%d after align.', kA);
    end
    if any(~cellfun(@(a,b) isequal(a,b), wins_in, winsL))
        warning('[LAP] Some windows were clipped to K=%d after DC/align.', kL);
    end

    % Bands per space, capped to K-1
    bandA = bandOpt; if isempty(bandA), bandA = max(3, round(0.06*kA)); end
    bandL = bandOpt; if isempty(bandL), bandL = max(3, round(0.06*kL)); end
    bandA = min(bandA, max(1, kA-1));
    bandL = min(bandL, max(1, kL-1));

    % ---- Compute windowed summaries (Adj / Lap) -----------------------------
    adj = compute_windows_global(VA_ref, VA_x, Lsym_ref, winsA, doSignAlign, SA_axis, modules, Wspatial, str, SArType, matchMode, bandA);
    lap = compute_windows_global(VL_ref, VL_x, Lsym_ref, winsL, doSignAlign, SA_axis, modules, Wspatial, str, SArType, matchMode, bandL, true);

    % ================== FIGURE 1: Windowed mean cosine ======================
    fig1 = myfig(['Spectral windows: ' str ' (mean cos (principal angles))'],[300 250 1020 460]);
    t = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    nexttile; 
        plot_window_bars(adj, winsA, cmap, labels, showLabels);
        title_if(showLabels, adj_label);
        grid on; set(gca,'XGrid','off','TickLength',[0 0]); box on; font(fntsz,usefont);
    nexttile; 
        plot_window_bars(lap, winsL, cmap, labels, showLabels);
        title_if(showLabels, lap_label);
        grid on; set(gca,'XGrid','off','TickLength',[0 0]); box on; font(fntsz,usefont);

    if showLabels
        ylabel(t,'Mean cosine','FontName',usefont,'FontSize',fntsz);
    end

    % ========== FIGURE 2+: Matched alignment per-k (with paging) ============
    K = size(adj.cos_k,1);                       
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
            grid on; set(gca,'XGrid','off'); box on; ylim([0 1]); font(fntsz,usefont);

        nexttile; 
            plot_perk_grouped(lap, cmap, clip_range(kShow, size(lap.cos_k,1)), showLabels); % safety if K differs
            title_if(showLabels, [lap_label ' — matched |cos| per k']);
            grid on; set(gca,'XGrid','off'); box on; ylim([0 1]); font(fntsz,usefont);

        if showLabels
            xlabel(t2,'Eigenvector index (k)'); 
            ylabel(t2,'|cos| matched','FontName',usefont,'FontSize',fntsz);
        end
    end

    % ---- Output -------------------------------------------------------------
    out = struct();
    out.adj = struct('Vref',VA_ref,'Vx',{VA_x},'eval_ref',evalA_ref,'eval_x',{evalA_x}, ...
                     'subspace_aff',{adj.subspace_aff}, 'match',{adj.match}, ...
                     'topology',{adj.topology}, 'perm',{adj.perm}, ...
                     'delta_names',{adj.delta_names});
    out.lap = struct('Vref',VL_ref,'Vx',{VL_x},'eval_ref',evalL_ref,'eval_x',{evalL_x}, ...
                     'subspace_aff',{lap.subspace_aff}, 'match',{lap.match}, ...
                     'topology',{lap.topology}, 'perm',{lap.perm}, ...
                     'delta_names',{lap.delta_names});
    out.figs = struct('subspace',fig1,'per_k',per_k_figs);
    out.labels=labels; out.cmap=cmap;
end

% ===================== CORE (global banded matching) ======================
function S = compute_windows_global(Vref, Vx_cell, Lref, wins, doSignAlign, SA_axis, modules, Wspatial, ~, SArType, matchMode, band, isLap)
    if nargin<15, isLap = false; end %#ok<NASGU>
    J = numel(Vx_cell);
    K = size(Vref,2);
    nW = numel(wins);

    % ensure in-range here too
    wins = clip_windows(wins, K);

    % Quick diagonal sign-align (optional)
    if doSignAlign
        for j=1:J, Vx_cell{j} = sign_align(Vref, Vx_cell{j}); end
    end

    % Global matching
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

    % Baseline diagonal cos
    cos_k = zeros(K,J);
    for j=1:J
        C = abs(Vref' * Vx_cell{j});
        cos_k(:,j) = diag(C);
    end

    % Window summaries + per-k topology (with Procrustes + sign-fix)
    subspace_aff = cell(1,nW);
    match = cell(1,nW);
    topology = cell(1,nW);
    delta_names = {}; % fill once

    for w=1:nW
        idx = wins{w};

        % Subspace affinity
        aff = zeros(1,J);
        for j=1:J
            vj = Vx_cell{j}(:, perm{j}(idx));
            aff(j) = subspace_affinity(Vref(:,idx), vj);
        end
        subspace_aff{w} = aff;

        % Matched + per-k metrics with Procrustes + sign-fix
        mcell = cell(1,J);
        topocell = cell(1,J);
        for j=1:J
            idxRef = idx(:)';                          
            idxX   = perm{j}(idxRef);                 
            idxX(idxX==0) = idxRef(idxX==0);

            cosMatched = zeros(1,numel(idxRef));
            T = repmat(struct('metrics',[],'delta',[]), 1, numel(idxRef));

            Vr = Vref(:,idxRef);
            Vx = Vx_cell{j}(:,idxX);
            [U,~,~] = svd(Vr'*Vx,'econ'); 
            R = U*U';           
            Vx_al = sign_fix(Vx*R, Vr, SA_axis, SArType);

            for t=1:numel(idxRef)
                r = idxRef(t); 
                vr = Vref(:,r);
                vx = Vx_al(:,t);

                cosMatched(t) = abs(vr' * vx);

                M = characterize_topology(vr, vx, Lref, SA_axis, modules, Wspatial, SArType);
                D = compute_deltas(M);
                T(t).metrics = M;
                T(t).delta   = D;

                if isempty(delta_names), delta_names = fieldnames(D); end
            end

            mcell{j} = struct('idxRef',idxRef,'idxX',idxX,'cosMatched',cosMatched);
            topocell{j} = T;
        end
        match{w} = mcell;
        topology{w} = topocell;
    end

    S = struct('subspace_aff',{subspace_aff}, 'match',{match}, 'topology',{topology}, ...
               'cos_k',cos_k, 'perm',{perm}, 'delta_names',{delta_names});
end

% ======================= Helpers (metrics & deltas) =======================
function metrics = characterize_topology(vr, vx, Lref, SA_axis, modules, Wspatial, SArType)
    metrics = struct();
    if ~isempty(SA_axis)
        metrics.SAr_ref = corr(vr, SA_axis, 'type',SArType,'rows','pairwise');
        metrics.SAr_x   = corr(vx, SA_axis, 'type',SArType,'rows','pairwise');
    else
        metrics.SAr_ref = NaN; metrics.SAr_x = NaN;
    end
    metrics.smooth_ref = max(eps, vr' * Lref * vr);
    metrics.smooth_x   = max(eps, vx' * Lref * vx);
    metrics.PR_ref = (sum(vr.^2)^2) / sum(vr.^4);
    metrics.PR_x   = (sum(vx.^2)^2) / sum(vx.^4);

    if ~isempty(modules)
        metrics.modR2_ref = module_R2(vr, modules);
        metrics.modR2_x   = module_R2(vx, modules);
    else
        metrics.modR2_ref = NaN; metrics.modR2_x = NaN;
    end

    if ~isempty(Wspatial)
        metrics.MoranI_ref = moransI(vr, Wspatial);
        metrics.MoranI_x   = moransI(vx, Wspatial);
    else
        metrics.MoranI_ref = NaN; metrics.MoranI_x = NaN;
    end
end

function D = compute_deltas(M)
    D = struct();
    if isfield(M,'SAr_x'),    D.delta_SAr    = M.SAr_x    - M.SAr_ref;    end
    if isfield(M,'smooth_x'), D.delta_smooth = M.smooth_x - M.smooth_ref; end
    if isfield(M,'PR_x'),     D.delta_PR     = M.PR_x     - M.PR_ref;     end
    if isfield(M,'modR2_x'),  D.delta_modR2  = M.modR2_x  - M.modR2_ref;  end
    if isfield(M,'MoranI_x'), D.delta_MoranI = M.MoranI_x - M.MoranI_ref; end
end

% =========================== Linear algebra utils =========================
function [Ahat, Lsym] = normalize_operators(A)
    A = symm(A);
    d = sum(A,2); d(d<=0) = eps;
    Dmh = spdiags(1./sqrt(d), 0, size(A,1), size(A,2));
    Ahat = Dmh * A * Dmh;
    Lsym = speye(size(A)) - Ahat;
end
function X = symm(X), X = (X+X')/2; end
function [V,D] = safe_eigs_laplacian(L, k)
    n = size(L,1); k = min(k, max(1,n-1));
    opts = struct('issym', true, 'isreal', true);
    try
        [V,D] = eigs((L+L')/2, k, 'smallestabs', opts);
    catch
        s = 1e-6 * max(1, full(trace(L))/n);
        [V,D] = eigs((L+L')/2, k, s, opts);
    end
end
function [V,e] = drop_dc_mode(V,e)
    N = size(V,1); c = ones(N,1)/sqrt(N);
    if abs(V(:,1)'*c) > 0.95
        V = V(:,2:end); e = e(2:end);
    end
end
function Vx = sign_align(Vref, Vx)
    C = Vref' * Vx; sgn = sign(diag(C)); sgn(sgn==0)=1;
    Vx = Vx .* sgn';
end
function V2 = sign_fix(Vx, Vr, SA, SArType)
    V2 = Vx;
    if isempty(SA), SA = zeros(size(Vx,1),1); end
    k = size(Vx,2);
    for i=1:k
        c1  = corr(Vx(:,i), Vr(:,i), 'type','pearson','rows','pairwise');
        c1b = corr(-Vx(:,i),Vr(:,i), 'type','pearson','rows','pairwise');
        c2  = corr(Vx(:,i), SA, 'type',SArType,'rows','pairwise');
        c2b = corr(-Vx(:,i),SA, 'type',SArType,'rows','pairwise');
        score_pos = nansum([abs(c1) + 0.5*abs(c2)]);
        score_neg = nansum([abs(c1b)+ 0.5*abs(c2b)]);
        if score_neg > score_pos, V2(:,i) = -Vx(:,i); end
    end
end
function aff = subspace_affinity(U, V)
    [~,S,~] = svd(U' * V, 'econ'); s = diag(S);
    aff = mean(min(max(s,0),1));
end
function perm = global_banded_match(Vref, Vx, band)
    K = size(Vref,2);
    S = abs(Vref' * Vx);  cost = 1 - S;
    [I,J] = ndgrid(1:K,1:K); %#ok<ASGLU>
    cost(abs(I-J) > band) = 1e6;
    if exist('matchpairs','file') == 2
        pairs = matchpairs(cost, 1e9);
        perm = zeros(K,1); perm(pairs(:,1)) = pairs(:,2);
    else
        assign = munkres(cost); perm = assign(:);
    end
    if any(perm==0)
        freeR = find(perm==0); usedC = perm(perm>0);
        avail = setdiff(1:K, usedC);
        for r = freeR(:)'
            [~,ix] = max(S(r,avail)); perm(r) = avail(ix); avail(ix) = [];
        end
    end
end

% ========================== Stats helpers =================================
function R2 = module_R2(v, mods)
    % One-way "module R^2": variance explained by module dummies
    v = v(:);
    u = unique(mods(:));
    X = zeros(numel(v), numel(u));
    for i = 1:numel(u)
        X(:,i) = double(mods == u(i));
    end
    b    = X \ v;
    vhat = X * b;

    sst = sum( (v - mean(v)).^2 );
    ssr = sum( (vhat - mean(v)).^2 );

    if sst <= eps
        R2 = 0;              % degenerate variance case
    else
        R2 = ssr / sst;      % explained / total
        R2 = min(1, max(0, R2));   % clamp to [0,1]
    end
end
function I = moransI(v, W)
    v = v(:); v0 = v - mean(v);
    num = v0' * W * v0; den = v0' * v0;
    I = num / max(den, eps);
end

% ========================== Plotting helpers ==============================
function plot_window_bars(S, wins, cmap, labels, showLabels)
    nW = numel(wins);
    J  = numel(S.subspace_aff{1});

    % Assemble matrix: rows = windows, cols = datasets
    M = zeros(nW, J);
    for w=1:nW, M(w,:) = S.subspace_aff{w}; end

    % Create grouped bars and keep handles
    B = bar(M, 'grouped', 'BarWidth', 0.75);   % B is 1xJ array of Bar objects

    % Color each dataset's bars (one object per dataset)
    if ~isempty(cmap) && numel(B) == J
        for j = 1:J
            B(j).FaceColor = 'flat';
            % nW-by-3 color rows (one row per window in this dataset)
            B(j).CData = repmat(cmap(j,:), nW, 1);
        end
    end

    % X tick labels = window ranges
    lbls = arrayfun(@(i) sprintf('%d-%d', wins{i}(1), wins{i}(end)), 1:nW, 'uni', 0);
    if showLabels
        set(gca, 'XTickLabel', lbls);
        if ~isempty(labels)
            legend(B, labels, 'Location','southoutside', 'Orientation','horizontal');
        end
    else
        set(gca,'XTickLabel',[]);
        legend off; title(''); xlabel(''); ylabel('');
    end
    ylim([0 1]);
end
function plot_perk_grouped(S, cmap, kShow, showLabels)
    if nargin<3 || isempty(kShow), kShow = []; end
    if nargin<4, showLabels = true; end

    K = size(S.cos_k,1);
    J = size(S.cos_k,2);

    % Start from baseline diagonal |cos|
    M = S.cos_k;

    % Overwrite with matched |cos| where available (per ref k)
    nW = numel(S.match);
    for w = 1:nW
        for j = 1:J
            m = S.match{w}{j};
            for t = 1:numel(m.idxRef)
                kr = m.idxRef(t);
                if kr >= 1 && kr <= size(M,1)
                    M(kr,j) = max(M(kr,j), m.cosMatched(t));
                end
            end
        end
    end

    if ~isempty(kShow)
        kShow = kShow(kShow >= 1 & kShow <= K);
        Mplot = M(kShow,:); xvals = kShow;
    else
        Mplot = M; xvals = 1:K;
    end

    % Grouped bar with handles
    B = bar(Mplot, 'grouped', 'BarWidth', 0.92);

    % Color each dataset object
    if ~isempty(cmap) && numel(B) == J
        for j = 1:J
            B(j).FaceColor = 'flat';
            B(j).CData = repmat(cmap(j,:), size(Mplot,1), 1);
        end
    end

    xlim([0 size(Mplot,1)+1]); ylim([0 1]);
    set(gca,'XTick',1:numel(xvals));
    if showLabels
        set(gca,'XTickLabel', xvals);
    else
        set(gca,'XTickLabel',[]);
        legend off; title(''); xlabel(''); ylabel('');
    end
end

% ========================== Tiny utils ====================================
function v = getOpt(s,f,d), if isfield(s,f), v=s.(f); else, v=d; end, end
function h = myfig(tit,pos), h=figure('Color','w','Position',pos,'Name',tit); end
function title_if(show, txt), if show, title(txt); else, title(''); end, end
function font(sz,fmly), set(gca,'FontSize',sz,'FontName',fmly); end
function arr = clip_range(arr, K), arr = arr(arr>=1 & arr<=K); end
function wins2 = clip_windows(wins, K)
    wins2 = wins;
    for i=1:numel(wins2)
        idx = wins2{i};
        idx = idx(idx >= 1 & idx <= K);
        if isempty(idx), idx = 1:min(3,K); end
        wins2{i} = idx;
    end
end
