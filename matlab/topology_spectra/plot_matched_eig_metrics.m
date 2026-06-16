function plot_matched_eig_metrics(out, operator, kSubset, focusWorst, topN, plotStyle, barMode, labelMode, groupMode)
% PLOT_MATCHED_EIG_METRICS
% Visualize topology metrics across matched eigenvectors for all datasets,
% and summarize the same metrics within each window (mean over selected ks).
%
% Inputs
%   out        : struct from spectral_alignment_windows_j()
%   operator   : 'adj' or 'lap'
%   kSubset    : optional vector of reference k indices to visualize
%   focusWorst : optional logical (default true if kSubset empty)
%   topN       : optional integer (default 2) number of worst ks per window
%   plotStyle  : optional 'heatmap' | 'bar' | 'both' (default: 'both')
%   barMode    : optional 'raw' | 'norm' (default: 'raw')
%   labelMode  : optional 'on' | 'off' (default: 'on')
%   groupMode  : optional 'window' (default) | 'metric' grouping of bar panels 
%
% Behavior
%   If kSubset provided -> use those ks.
%   Else if focusWorst (default) -> identify, for each window, the topN worst-aligned
%   eigs (lowest matched |cos|) pooled across datasets, and union them.
%
% Displays
%   If 'heatmap' or 'both':
%     Figure 1: heatmaps of Δ(X−Ref) per-k (rows=datasets, cols=k indices).
%     Figure 2: heatmaps of window-mean Δ(X−Ref) (rows=datasets, cols=windows).
%   If 'bar' or 'both':
%     Figure 3: per-window grouped **bar** summaries (x=metrics, grouped by dataset).
%
% 2025 Mark C Nelson (MNI)
% -------------------------------------------------------------------------

    if nargin < 6 || isempty(plotStyle), plotStyle = 'both';    end
    if nargin < 7 || isempty(barMode),   barMode   = 'raw';     end
    if nargin < 8 || isempty(labelMode), labelMode = 'on';      end
    if nargin < 9 || isempty(groupMode), groupMode = 'window';  end         % governs how bar plot panels are grouped
    % plotStyle = lower(plotStyle);
    % barMode   = lower(barMode);
    % labelMode = lower(labelMode);

    cm = colormap(smartcmaps('bentcoolwarm'));
    fntsz=18; fntname='Helvetica';

    R = out.(operator);                 % .match, .topology, (and usually .perm)
    J = numel(R.Vx);
    wins = R.match;                     % cell{w}{j} with fields .idxRef,...
    topo = R.topology;                  % cell{w}{j}(t).metrics
    nW = numel(wins);

    % dataset labels & colors
    if isfield(out,'labels') && numel(out.labels)==J
        dsLabels = out.labels(:).';
    else
        dsLabels = arrayfun(@(j) sprintf('D%d',j), 1:J, 'UniformOutput', false);
    end
    if isfield(out,'cmap') && size(out.cmap,1)==J && size(out.cmap,2)==3
        dsColors = out.cmap;
    else
        dsColors = lines(J);
    end

    % ---- Determine k set
    if nargin < 3 || isempty(kSubset)
        if nargin < 4 || isempty(focusWorst), focusWorst = true; end
        if nargin < 5 || isempty(topN),       topN = 2;         end
        if focusWorst
            kSel = [];
            for w=1:nW
                Klist = []; Coslist = [];
                for j=1:J
                    m = wins{w}{j};
                    for t=1:numel(m.idxRef)
                        Klist(end+1)   = m.idxRef(t); %#ok<AGROW>
                        Coslist(end+1) = m.cosMatched(t); %#ok<AGROW>
                    end
                end
                [uk,~,ic] = unique(Klist);
                mincos = accumarray(ic(:), Coslist(:), [], @min);
                [~,ord] = sort(mincos, 'ascend');
                take = uk(ord(1:min(topN, numel(uk))));
                kSel = union(kSel, take);
            end
            kSubset = sort(unique(kSel));
        else
            error('Provide kSubset or enable focusWorst.');
        end
    end

    % ---- Collect metrics into Δ matrices: J x Ksel
    Ksel = numel(kSubset);
    Mats = struct();  % each field -> JxKsel matrix

    fieldsIn = {'SA_r','smooth','PR','MoranI','modR2'};
    for j=1:J
        for kk=1:Ksel
            k = kSubset(kk);
            % find which window (for dataset j) contains k
            wHit = 0; tHit = 0;
            for w=1:nW
                m = wins{w}{j};
                t = find(m.idxRef==k, 1);
                if ~isempty(t), wHit = w; tHit = t; break; end
            end
            if wHit==0, continue; end
            metr = topo{wHit}{j}(tHit).metrics;

            % fill Mats.*(j,kk)
            if isfield(metr,'SAr_ref') && isfield(metr,'SAr_x')
                if ~isfield(Mats,'SA_r'), Mats.SA_r = nan(J,Ksel); end
                Mats.SA_r(j,kk) = metr.SAr_x - metr.SAr_ref;
            end
            if isfield(metr,'smooth_ref') && isfield(metr,'smooth_x')
                if ~isfield(Mats,'smooth'), Mats.smooth = nan(J,Ksel); end
                Mats.smooth(j,kk) = metr.smooth_x - metr.smooth_ref;
            end
            if isfield(metr,'PR_ref') && isfield(metr,'PR_x')
                if ~isfield(Mats,'PR'), Mats.PR = nan(J,Ksel); end
                Mats.PR(j,kk) = metr.PR_x - metr.PR_ref;
            end
            if isfield(metr,'MoranI_ref') && isfield(metr,'MoranI_x')
                if ~isfield(Mats,'MoranI'), Mats.MoranI = nan(J,Ksel); end
                Mats.MoranI(j,kk) = metr.MoranI_x - metr.MoranI_ref;
            end
            if isfield(metr,'modR2_ref') && isfield(metr,'modR2_x')
                if ~isfield(Mats,'modR2'), Mats.modR2 = nan(J,Ksel); end
                Mats.modR2(j,kk) = metr.modR2_x - metr.modR2_ref;
            end
        end
    end

    % Which metrics actually exist?
    names = fieldnames(Mats);
    if isempty(names)
        warning('No metrics available to plot. Did you provide SA_axis/modules/Wspatial?');
        return;
    end

    % ---- Figure 1: per-k heatmaps
    if any(strcmpi(plotStyle, {'both','heatmap'}))
        fig1 = myfig(sprintf('Matched-eig metrics (%s) — Δ(X−Ref): %s', upper(operator), mat2str(kSubset)), [240 120 1120 720]);
        T1 = tiledlayout(numel(names),1,'Padding','compact','TileSpacing','compact');
        for i=1:numel(names)
            M = Mats.(names{i});                     % J x Ksel
            nexttile;
            imagesc(1:Ksel, 1:J, M);                 % plot in index space
            axis tight; colorbar;
            cax = clim; CMAX = max(abs(cax)); clim([-CMAX CMAX]);
            title([pretty_name(names{i}) '  (X − Ref)']);
            set(gca,'YDir','normal'); font(fntsz,fntname);
            set(gca,'YTick',1:J,'YTickLabel',dsLabels);
            [xt, xl] = make_sparse_ticks(1:Ksel, kSubset, 10);
            set(gca,'XTick',xt,'XTickLabel',xl);
        end
        title(T1, sprintf('Operator: %s — Δ metrics across matched eigenvectors', upper(operator)));
        colormap(smartcmaps('bentcoolwarm'));
    end

    % ================== WINDOW SUMMARIES (heatmaps) ==================
    % window index lists from reference indices (use dataset 1 as template)
    winKs = cell(1,nW);
    for w=1:nW
        winKs{w} = wins{w}{1}.idxRef(:)';   % global ref k in this window
    end
    winLbls = arrayfun(@(w) sprintf('k %d–%d', min(winKs{w}), max(winKs{w})), 1:nW, 'UniformOutput', false);

    % Compute J x nW means for each metric over selected ks within the window
    MatsWin = struct();  % each field -> J x nW
    for i=1:numel(names)
        M = Mats.(names{i});              % J x Ksel
        MW = nan(J, nW);
        for w=1:nW
            ks = intersect(kSubset, winKs{w});
            if isempty(ks)
                MW(:,w) = nan;
            else
                [~, loc] = ismember(ks, kSubset);
                MW(:,w) = mean(M(:, loc), 2, 'omitnan');
            end
        end
        MatsWin.(names{i}) = MW;
    end

    if any(strcmpi(plotStyle, {'both','heatmap'}))
        fig2 = myfig(sprintf('Matched-eig metrics (%s) — WINDOW MEANS of Δ(X−Ref)', upper(operator)), [260 140 1120 720]);
        T2 = tiledlayout(numel(names),1,'Padding','compact','TileSpacing','compact');
        for i=1:numel(names)
            MW = MatsWin.(names{i});                 % J x nW
            nexttile;
            imagesc(1:nW, 1:J, MW);
            axis tight; colorbar;
            cax = clim; CMAX = max(abs(cax)); clim([-CMAX CMAX]);
            title([pretty_name(names{i}) '  — window mean Δ( X − Ref )']);
            set(gca,'YDir','normal'); font(fntsz,fntname);
            set(gca,'XTick',1:nW,'XTickLabel',winLbls,'YTick',1:J,'YTickLabel',dsLabels);
        end
        title(T2, sprintf('Operator: %s — windowed means over displayed k''s', upper(operator)));
        colormap(smartcmaps('bentcoolwarm'));
    end

    % ================== WINDOW SUMMARIES (BAR PLOTS) ==================
    if any(strcmpi(plotStyle, {'both','bar'}))

        % Assemble data cube: [J datasets] x [numMetrics] x [nW windows]
        numM = numel(names);
        Cube = nan(J, numM, nW);
        for i=1:numM
            Cube(:, i, :) = MatsWin.(names{i}); % J x 1 x nW -> J x nW along 3rd dim
        end

        % If normalized mode, scale each metric independently by its global max |Δ|
        scaleText = '';
        if strcmpi(barMode,'norm')
            for i=1:numM
                X = squeeze(Cube(:,i,:));           % J x nW
                s = max(abs(X(:)), [], 'omitnan');
                if s>0, Cube(:,i,:) = Cube(:,i,:) ./ s; end
            end
            scaleText = ' (relative)';
        end

        % Decide how bars should be grouped
        if strcmpi(groupMode,'metric')
            dim_in   = nW;     labels_in   = winLbls;   
            dim_out  = numM;   labels_out  = names;
        else                                                                % else default to by window organization
            dim_out  = nW;     labels_out  = winLbls;   
            dim_in   = numM;   labels_in   = names;
        end

        % Build a single 1×nW figure; each tile = a window, x-groups = metrics, bars = datasets
        fig3 = myfig(sprintf('Matched-eig metrics (%s) — WINDOW MEANS (bars, %s)', upper(operator), barMode), [260 160 1200 520]);
        T3 = tiledlayout(1, dim_out, 'Padding','compact', 'TileSpacing','compact');

        for w=1:dim_out
            nexttile; hold on;
            % Either a J x numM or J x nW matrix (bar expects rows = groups)
            if strcmpi(groupMode,'metric')
                MWmat = squeeze(Cube(:,w,:)).';                             % J x nW                    
            else
                MWmat = squeeze(Cube(:,:,w)).';                             % J x numM
            end
            B = bar(MWmat, 'grouped', 'BarWidth', 0.75);
            if strcmpi(barMode,'norm'); ylim([-1 1]); end
            for j=1:J
                B(j).FaceColor = 'flat';
                B(j).CData = repmat(dsColors(j,:), dim_in, 1);
            end
            yline(0,'k-','LineWidth',0.75);
            set(gca,'XTick',1:dim_in,'XTickLabel','','TickLength',[0 0]);
            grid on; set(gca,'XGrid','off'); box on; font(fntsz,fntname);
            for ll=1:1:dim_in-1; line([ll+.5 ll+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end
            if strcmpi(labelMode,'on')
                title(sprintf('%s', labels_out{w}));
                set(gca,'XTickLabel',cellfun(@pretty_name, names, 'UniformOutput', false)); xtickangle(30);
                if w==dim_out; legend(B, dsLabels, 'Location','northoutside', 'Orientation','horizontal'); end
                if w==1;  ylabel(sprintf('Δ(X−Ref)%s', scaleText)); end
            else
                if w~=1; set(gca,'YTickLabel',''); end
            end
            hold off;
        end
        % title(T3, sprintf('Operator: %s — window means (bars: metrics grouped; colors = datasets)', upper(operator)));
    end
end

% ---------- helpers ----------
function s = pretty_name(f)
    switch lower(f)
        case 'sa_r',   s = 'SA_r';
        case 'smooth', s = 'smooth';
        case 'pr',     s = 'PR';
        case 'morani', s = 'MoranI';
        case 'modr2',  s = 'modR2';
        otherwise,     s = f;
    end
end

function [xticks, xlabels] = make_sparse_ticks(xidx, xlabels_full, maxN)
    if nargin < 3 || isempty(maxN), maxN = 10; end
    K = numel(xidx);
    if K <= maxN
        xticks = xidx;
        xlabels = arrayfun(@num2str, xlabels_full, 'UniformOutput', false);
    else
        pos = unique(round(linspace(1, K, maxN)));
        xticks = xidx(pos);
        xlabels = arrayfun(@(v) num2str(v), xlabels_full(pos), 'UniformOutput', false);
    end
end
