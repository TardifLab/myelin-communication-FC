function figs = delta_topology_boxcharts(out, which_space, metric_list, labels, cmap, labelMode)
% DELTA_TOPOLOGY_BOXCHARTS (normalized to [-1,1] with y=0 baseline)
% Window-wise boxcharts of per-k delta topology metrics (x - ref) per dataset.
% - All metrics are robust-normalized to [-1, 1] using the 99th percentile
%   of absolute values pooled across all windows/datasets for that metric,
%   then clipped to [-1,1].
% - Common ylim = [-1, 1] on every panel; thin black baseline at y=0.
%
% Usage:
%   figs = delta_topology_boxcharts(out,'lap', ...
%          {'delta_SAr','delta_smooth','delta_PR','delta_modR2','delta_MoranI'}, ...
%          labels, cmap, 'on');
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
% -------------------------------------------------------------------------

if nargin < 6 || isempty(labelMode), labelMode = 'on';   end
labelMode = lower(labelMode);
fntsz=20; fntname='Calibri';


assert(isfield(out,which_space),'which_space must be ''adj'' or ''lap''.');
S = out.(which_space);

if nargin < 3 || isempty(metric_list), metric_list = S.delta_names; end
if nargin < 4 || isempty(labels),      labels = out.labels;          end
if nargin < 5 || isempty(cmap),        cmap   = out.cmap;            end

nW = numel(S.topology);
J  = numel(S.topology{1});
figs = gobjects(numel(metric_list),1);

for m = 1:numel(metric_list)
    met = metric_list{m};

    % Gather raw deltas {w}{j}
    Delta = cell(nW,J);
    for w = 1:nW
        for j = 1:J
            T = S.topology{w}{j};
            vals = nan(numel(T),1);
            for t = 1:numel(T)
                if isfield(T(t),'delta') && isfield(T(t).delta, met)
                    vals(t) = T(t).delta.(met);
                end
            end
            Delta{w,j} = vals(~isnan(vals));
        end
    end

    % ---- Robust, sign-preserving normalization to [-1,1] ----
    allY = vertcat(Delta{:});               % pool across windows & datasets
    allY = allY(~isnan(allY));
    if isempty(allY)
        scale = 1;  % nothing to scale; avoid divide-by-zero
    else
        % Use 99th percentile of |values| for robust scaling, keep sign
        scale = prctile(abs(allY), 99);
        if ~isfinite(scale) || scale <= 0, scale = max(1, max(abs(allY))); end
    end
    for w = 1:nW
        for j = 1:J
            y = Delta{w,j};
            if isempty(y), continue; end
            y = y ./ scale;                 % scale
            y = max(-1, min(1, y));         % clip
            Delta{w,j} = y;
        end
    end

    % ---- Figure (one row; columns = windows) ----
    figs(m) = myfig(sprintf('%s — %s (normalized to [-1,1])', upper(which_space), met),[220 160 1180 540]);
    tl = tiledlayout(1,nW,'Padding','compact','TileSpacing','compact');

    for w = 1:nW
        nexttile; hold on

        % Plot one box per dataset at x=j
        for j = 1:J
            y = Delta{w,j};
            if isempty(y), continue; end
            x = j * ones(numel(y),1);

            boxchart(x, y, ...
                'BoxFaceColor', cmap(j,:), ...
                'BoxFaceAlpha', 0.7, ...
                'BoxEdgeColor', cmap(j,:), ...
                'WhiskerLineColor', cmap(j,:), ...
                'MarkerStyle','o', ...
                'LineWidth', 2, ...
                'BoxMedianLineColor','k');                   % thickens box edges, whiskers, median
        end

        % Axes cosmetics + common limits and baseline
        xlim([0.5, J+0.5]); grid on; set(gca,'XGrid','off'); box on; ylim([-1 1]);  % common y-limits
        yline(0,'k-','LineWidth',0.75);                                             % baseline at zero
        font(fntsz,fntname);

      % Tile labels
        if strcmp(labelMode,'on')
            title(sprintf('win %d–%d', S.match{w}{1}.idxRef(1), S.match{w}{1}.idxRef(end)));
            set(gca,'XTick',1:J,'XTickLabel',labels,'XTickLabelRotation',20);
          % Do these only on 1st tile
            if w == 1
                ylabel(sprintf('%s', strrep(met,'_','\_')));
                ph = gobjects(J,1);
                for j = 1:J
                    ph(j) = patch(NaN,NaN,cmap(j,:), 'EdgeColor',cmap(j,:), 'FaceAlpha',0.7);
                end
                legend(ph, labels, 'Location','southoutside','Orientation','horizontal'); % (color patches), once on first tile
            end
        else
            set(gca,'XTick',1:J,'XTickLabel','');
            if w~=1; set(gca,'YTickLabel',''); end % turn of yticklabels for all but 1st tile
        end
    end
end
%--------------------------------------------------------------------------
end
