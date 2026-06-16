function out = analyze_spectral_alignment_windows(A_ref,A_x,str,opts)
% Analyze spectral alignment of target network A_x vs reference A_ref
% across adjacency & Laplacian operators using windowed subspaces, with
% optional kernel-weighted overlays and eigenvector topology characterization.
%
% Inputs
%   A_ref : NxN weighted adjacency (reference; e.g., Caliber)
%   A_x   : NxN weighted adjacency (target; e.g., Delay or MTsat)
%   str   : char vector for plot labels
%   opts (struct with fields; all optional)
%     .kmax            (default 50)
%     .normalize       (true/false) -> use Ahat = D^(-1/2) A D^(-1/2), Lsym = I - Ahat
%     .wins            (cell array of index ranges), e.g. {2:5, 5:10, 10:15, 15:25, 25:50}
%     .doSignAlign     (default true) sign-align matched eigenvectors
%     .kernel.doAdj    (default false) overlay communicability weights for beta
%     .kernel.doLap    (default false) overlay heat-kernel weights for tau
%     .kernel.beta     (default 1) communicability scale (rescaled internally)
%     .kernel.tau      (default 1) heat-kernel scale (rescaled internally)
%     .coords          (Nx3) node coords for spatial map (optional)
%     .SA_axis         (Nx1) sensory-association gradient (optional)
%     .modules         (Nx1 int) module labels (optional)
%     .Wspatial        (NxN) spatial weight matrix for Moran's I (row-normalized recommended)
%     .plotMaps        (true/false) plot matched eigenvector maps per window (default false)
%
% Output (struct)
%   .adj, .lap : per-operator results with fields:
%       .Vref, .Vx, .eval_ref, .eval_x     (eigenvectors/values)
%       .subspace_aff(win)                 (mean cos principal angles per window)
%       .match(win).idxRef, .idxX          (matched indices within window)
%       .match(win).meanCorr               (mean matched correlation in window)
%       .match(win).topology(win, k).metrics (if requested): SA r, smoothness, PR, MoranI, module R2
%   .figs : handles of summary figures
%
% Requires: Statistics & Machine Learning Toolbox (matchpairs).
% If unavailable, replace matchpairs() with a Hungarian (munkres) implementation.
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------


% ---------------- Defaults ----------------
if nargin < 3, opts = struct(); end
kmax        = getOpt(opts, 'kmax', 50);
normalize   = getOpt(opts, 'normalize', true);
wins        = getOpt(opts, 'wins', {2:5, 5:10, 10:15, 15:25, 25:50});
doSignAlign = getOpt(opts, 'doSignAlign', true);
coords      = getOpt(opts, 'coords', []);
SA_axis     = getOpt(opts, 'SA_axis', []);
modules     = getOpt(opts, 'modules', []);
Wspatial    = getOpt(opts, 'Wspatial', []);
plotMaps    = getOpt(opts, 'plotMaps', false);

kernel.doAdj = false; kernel.doLap = false; kernel.beta = 1; kernel.tau = 1;
if isfield(opts,'kernel'), kernel = mergeStruct(kernel, opts.kernel); end

% ------------- Build operators -------------
if normalize
    [Ahat_ref, Lsym_ref] = normalize_operators(A_ref);
    [Ahat_x,   Lsym_x]   = normalize_operators(A_x);
    adj_label = 'Adj (norm)'; lap_label = 'Lap (norm)';
else
    Ahat_ref = symm(A_ref); Ahat_x = symm(A_x);
    Lsym_ref = diag(sum(A_ref,2)) - Ahat_ref;
    Lsym_x   = diag(sum(A_x,2))   - Ahat_x;
    adj_label = 'Adj'; lap_label = 'Lap';
end

% ------------- Eigendecompositions ----------
[VA_ref, DA_ref] = eigs(Ahat_ref, kmax, 'largestabs');
[VA_x,   DA_x]   = eigs(Ahat_x,   kmax, 'largestabs');
[VL_ref, DL_ref] = eigs(Lsym_ref, kmax, 'smallestabs');
[VL_x,   DL_x]   = eigs(Lsym_x,   kmax, 'smallestabs');

[VA_ref, VA_x] = norm_and_sign(VA_ref, VA_x, doSignAlign);
[VL_ref, VL_x] = norm_and_sign(VL_ref, VL_x, doSignAlign);

evalA_ref = sort(diag(DA_ref)); evalA_x = sort(diag(DA_x));
evalL_ref = sort(diag(DL_ref)); evalL_x = sort(diag(DL_x));

% ------------- Windowed subspace alignment & matching -------------
adj = compute_windows(VA_ref, VA_x, wins);
lap = compute_windows(VL_ref, VL_x, wins);

% ------------- Kernel-weighted contributions (optional) -------------
if kernel.doAdj
    % scale beta relative to spectral radius so beta*rho is constant across networks
    rho = max(abs(diag(DA_ref)));
    betaEff = kernel.beta / max(rho, eps);
    wA = mode_weights(diag(DA_ref), 'adj', betaEff);
else
    wA = [];
end

if kernel.doLap
    lamMax = max(diag(DL_ref));
    tauEff = kernel.tau / max(lamMax, eps);
    wL = mode_weights(diag(DL_ref), 'lap', tauEff);
else
    wL = [];
end

% ------------- Topology characterization of matched modes -------------
doTopo = ~isempty(coords) || ~isempty(SA_axis) || ~isempty(modules) || ~isempty(Wspatial);
if doTopo
    adj = characterize_windows_topology(adj, VA_ref, VA_x, Lsym_ref, coords, SA_axis, modules, Wspatial, plotMaps, 'Adjacency', str);
    lap = characterize_windows_topology(lap, VL_ref, VL_x, Lsym_ref, coords, SA_axis, modules, Wspatial, plotMaps, 'Laplacian', str);
end

% ------------- Plots: windowed subspace + (optional) kernel overlays ---
fig1 = myfig(['Comparing subspace: ' str],[200 200 1100 480]);
subplot(1,2,1);
bar(cellfun(@(c)c, adj.subspace_aff), 0.8); ylim([0 1]); grid on;
set(gca,'XTick',1:numel(wins),'XTickLabel',winlabels(wins)); xtickangle(30);
ylabel('Mean cos(principal angles)'); title(['Subspace alignment - ' adj_label]); font(20,'Calibri');

subplot(1,2,2);
bar(cellfun(@(c)c, lap.subspace_aff), 0.8); ylim([0 1]); grid on;
set(gca,'XTick',1:numel(wins),'XTickLabel',winlabels(wins)); xtickangle(30);
ylabel('Mean cos(principal angles)'); title(['Subspace alignment - ' lap_label]); font(20,'Calibri');

% Optional kernel overlays (add small inset lines weighted by wA/wL)
if kernel.doAdj || kernel.doLap
    fig2 = myfig(['Alignment & kernel contribution: ' str],[200 200 1100 480]); 
    if kernel.doAdj
        subplot(1,2,1);
        plot(1:kmax, abs(diag(VA_ref'*VA_x)), 'k-', 'LineWidth', 2); hold on;
        yyaxis right; plot(1:kmax, wA( order_by_ref(VA_ref, Ahat_ref) ), '.', 'MarkerSize', 20);
        ylabel('Mode weight'); yyaxis left; ylim([0 1]);
        xlabel('k'); ylabel('Cosine'); title(['Adj alignment & kernel weights (\beta=' num2str(kernel.beta) ')']); grid on;
        legend({'cos(v_{ref,k}, v_{x,k})','kernel weight'}, 'Location','best'); hold off; font(16,'Calibri');
    end
    if kernel.doLap
        subplot(1,2,2);
        plot(1:kmax, abs(diag(VL_ref'*VL_x)), 'k-', 'LineWidth', 2); hold on;
        yyaxis right; plot(1:kmax, wL( order_by_ref(VL_ref, Lsym_ref) ), '.', 'MarkerSize', 20);
        ylabel('Mode weight'); yyaxis left; ylim([0 1]);
        xlabel('k'); ylabel('Cosine'); title(['Lap alignment & kernel weights (\tau=' num2str(kernel.tau) ')']); grid on;
        legend({'cos(v_{ref,k}, v_{x,k})','kernel weight'}, 'Location','best'); hold off; font(16,'Calibri');
    end
else
    fig2 = [];
end

% ------------- Pack outputs -------------
out = struct();
out.adj = struct('Vref',VA_ref,'Vx',VA_x,'eval_ref',evalA_ref,'eval_x',evalA_x, ...
                 'subspace_aff',{adj.subspace_aff}, 'match',{adj.match});
out.lap = struct('Vref',VL_ref,'Vx',VL_x,'eval_ref',evalL_ref,'eval_x',evalL_x, ...
                 'subspace_aff',{lap.subspace_aff}, 'match',{lap.match});
out.figs = struct('subspace',fig1,'kernel',fig2);

end

% ======================= Helpers =======================

function X = symm(X), X = (X+X')/2; end
function v = getOpt(s, f, d), if isfield(s,f), v=s.(f); else, v=d; end
end

function S = mergeStruct(S, T)
f = fieldnames(T);
for i=1:numel(f), S.(f{i}) = T.(f{i}); end
end

function [Ahat, Lsym] = normalize_operators(A)
A = (A + A')/2;
s = sum(A,2);
invS = 1 ./ max(sqrt(s), eps);
Dmh = diag(invS);
Ahat = Dmh * A * Dmh;       % normalized adjacency in [-1,1]
Lsym = eye(size(A)) - Ahat; % normalized Laplacian in [0,2]
end

function [Vref, Vx] = norm_and_sign(Vref, Vx, doSign)
Vref = bsxfun(@rdivide, Vref, sqrt(sum(Vref.^2,1)));
Vx   = bsxfun(@rdivide, Vx,   sqrt(sum(Vx.^2,1)));
if doSign
    for i=1:size(Vref,2)
        if corr(Vref(:,i), Vx(:,i)) < 0, Vx(:,i) = -Vx(:,i); end
    end
end
end

function labels = winlabels(wins)
labels = cellfun(@(w) sprintf('%d-%d', w(1), w(end)), wins, 'UniformOutput', false);
end

function res = compute_windows(Vref, Vx, wins)
res.subspace_aff = cell(numel(wins),1);
res.match        = cell(numel(wins),1);
for w = 1:numel(wins)
    idx = wins{w};
    % Subspace affinity via principal angles
    [~,S,~] = svd(Vref(:,idx)' * Vx(:,idx), 'econ');
    res.subspace_aff{w} = mean(diag(S)); % mean cos(theta)
    % One-to-one matching within window
    C = abs(Vref(:,idx)' * Vx(:,idx));   % pairwise cosines
    % Try matchpairs (Stats TB); otherwise greedy matching
    try
        pairs = matchpairs(1 - C, 1e9); % minimize (1-C) => maximize C
        meanCorr = mean(arrayfun(@(r) C(r, pairs(r,2)), 1:size(C,1)));
        assign = pairs(:,2)';            % column indices in Vx for each idx row
    catch
        % greedy fallback
        assign = zeros(1, numel(idx));
        Ctmp = C;
        for k = 1:numel(idx)
            [m,ind] = max(Ctmp(:));
            [r,c] = ind2sub(size(Ctmp), ind);
            assign(r) = c; Ctmp(r,:) = -inf; Ctmp(:,c) = -inf;
        end
        meanCorr = mean(C(sub2ind(size(C), 1:numel(idx), assign)));
    end
    res.match{w} = struct('idxRef', idx(:)', 'idxX', idx(assign), 'meanCorr', meanCorr);
end
end

function ord = order_by_ref(Vref, ~)
% placeholder: identity ordering by index; customize if needed
ord = 1:size(Vref,2);
end

function w = mode_weights(eigvals_diag, type, param)
lam = eigvals_diag;  % vector
switch lower(type)
    case 'adj'
        w = exp(param * lam);
    case 'lap'
        w = exp(-param * lam);
    otherwise
        error('type must be ''adj'' or ''lap''');
end
w = w / sum(w);
end

function res = characterize_windows_topology(res, Vref, Vx, Lref, coords, SA, modules, Wspatial, doMaps, tag, str)
% For each window, compute topology metrics for matched eigenvectors.
% Metrics per matched mode: SA corr, smoothness v^T L v, participation ratio,
% module R^2 (from one-hot regression), Moran's I (if Wspatial).
for w = 1:numel(res.match)
    m = res.match{w};
    K = numel(m.idxRef);
    topo = struct([]);
    for k = 1:K
        vr = Vref(:, m.idxRef(k));
        vx = Vx(:,   m.idxX(k));
        metrics = struct();

        % SA-axis correlation (if provided)
        if ~isempty(SA)
            metrics.SA_r_ref = corr(vr, SA, 'type','Spearman','rows','pairwise');
            metrics.SA_r_x   = corr(vx, SA, 'type','Spearman','rows','pairwise');
        end

        % Smoothness wrt reference Laplacian
        metrics.smooth_ref = vr' * Lref * vr;
        metrics.smooth_x   = vx' * Lref * vx;

        % Participation Ratio (localization)
        metrics.PR_ref = (sum(vr.^2)^2) / sum(vr.^4);
        metrics.PR_x   = (sum(vx.^2)^2) / sum(vx.^4);

        % Module expressivity (R^2 of GLM with module dummies)
        if ~isempty(modules)
            metrics.modR2_ref = module_R2(vr, modules);
            metrics.modR2_x   = module_R2(vx, modules);
        end

        % Moran's I (spatial autocorr)
        if ~isempty(Wspatial)
            metrics.MoranI_ref = moransI(vr, Wspatial);
            metrics.MoranI_x   = moransI(vx, Wspatial);
        end

        topo(k).metrics = metrics;

        % Optional simple spatial maps
        if doMaps && ~isempty(coords)
            myfig(['Spatial maps of eigenvectors: ' str],[200 200 920 720]);

            % --- Custom red-white-blue colormap centered on 0
            % % % n = 256;
            % % % half = floor(n/2);
            % % % blue2white = [linspace(0,1,half)'  linspace(0,1,half)'  ones(half,1)];
            % % % white2red  = [ones(n-half,1)       linspace(1,0,n-half)' linspace(1,0,n-half)'];
            % % % cmap_rwb = [blue2white; white2red];
            % % % cmap_rwb(half-2:half+2,:) = 1;                     % narrow white seam near 0
            cmap_rwb=smartcmaps('bentcoolwarm');
            
            % Compute color limits symmetric around 0
            CLIM = max(abs([vr(:); vx(:)]));
            CLIM = [-CLIM CLIM];

            % VIEW 1: lateral
            subplot(2,2,1);
            scatter3(coords(:,1),coords(:,2),coords(:,3),24,vr,'filled');
            axis equal off; view([90 0]);  % adjust as desired (e.g. sagittal)
            title(sprintf('%s win %d | ref k=%d', tag(1:3), w, m.idxRef(k)));
            colormap(cmap_rwb); colorbar; clim(CLIM); font(18,'Calibri');
            
            subplot(2,2,2);
            scatter3(coords(:,1),coords(:,2),coords(:,3),24,vx,'filled');
            axis equal off; view([90 0]);
            title(sprintf('%s win %d | x k=%d', tag(1:3), w, m.idxX(k)));
            colormap(cmap_rwb); colorbar; clim(CLIM); font(18,'Calibri');
            
            % VIEW 2: axial
            subplot(2,2,3);
            scatter3(coords(:,1),coords(:,2),coords(:,3),24,vr,'filled');
            axis equal off; view([0 90]);  % top-down view
            title(sprintf('%s win %d | ref k=%d', tag(1:3), w, m.idxRef(k)));
            colormap(cmap_rwb); colorbar; clim(CLIM); font(18,'Calibri');

            subplot(2,2,4);
            scatter3(coords(:,1),coords(:,2),coords(:,3),24,vx,'filled');
            axis equal off; view([0 90]);
            title(sprintf('%s win %d | x k=%d', tag(1:3), w, m.idxX(k)));
            colormap(cmap_rwb); colorbar; clim(CLIM); font(18,'Calibri');

        end
    end
    res.match{w}.topology = topo;
end
end

function R2 = module_R2(v, modules)
% Regress v on module dummies, return R^2 (variance explained by modules)
mods = modules(:); u = unique(mods(:));
X = zeros(numel(v), numel(u));
for i=1:numel(u), X(:,i) = double(mods==u(i)); end
b = X \ v; vhat = X*b;
ssr = sum((vhat-mean(v)).^2); sst = sum((v-mean(v)).^2);
R2 = max(0, min(1, ssr / max(sst, eps)));
end

function I = moransI(v, W)
% Moran's I with row-stochastic W (or any symmetric W)
v = v(:); v0 = v - mean(v);
num = v0' * W * v0;
den = v0' * v0;
I = num / max(den, eps);
end
