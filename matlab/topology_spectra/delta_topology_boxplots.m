function figs = delta_topology_boxplots(out, which_space, metric_list, labels, cmap)
% DELTA_TOPOLOGY_BOXPLOTS
% Make window-wise boxplots (median line; whiskers to Q1/Q3) of per-k
% delta topology metrics (x - ref) across all k in each window, per dataset.
%
% Inputs:
%   out          : struct from spectral_alignment_windows_j
%   which_space  : 'adj' or 'lap'
%   metric_list  : cellstr of delta metric names to plot, e.g.,
%                  {'delta_SAr','delta_smooth','delta_PR','delta_modR2','delta_MoranI'}
%   labels       : 1xJ cellstr dataset labels
%   cmap         : Jx3 colors
%
% Output:
%   figs         : array of figure handles (one per metric)
%
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
% -------------------------------------------------------------------------

    assert(isfield(out,which_space), 'which_space must be ''adj'' or ''lap''.');
    S = out.(which_space);

    if nargin<3 || isempty(metric_list)
        metric_list = S.delta_names; % whatever was detected in compute_deltas
    end
    if nargin<4 || isempty(labels), labels = out.labels; end
    if nargin<5 || isempty(cmap),    cmap    = out.cmap;  end

    nW = numel(S.topology);
    J  = numel(S.topology{1});
    figs = gobjects(numel(metric_list),1);

    % Gather deltas per metric/window/dataset
    for m = 1:numel(metric_list)
        met = metric_list{m};

        % Build a cell {w}{j} -> vector of deltas across k in window w
        Delta = cell(nW,J);
        for w=1:nW
            for j=1:J
                T = S.topology{w}{j};
                vec = nan(1, numel(T));
                for t=1:numel(T)
                    if isfield(T(t),'delta') && isfield(T(t).delta, met)
                        vec(t) = T(t).delta.(met);
                    end
                end
                Delta{w,j} = vec(:);
            end
        end

        % Plot: one figure per metric; panels = windows; grouped by dataset
        figs(m) = figure('Color','w','Position',[200 180 1120 520], ...
                         'Name',sprintf('%s — %s', upper(which_space), met));
        t = tiledlayout(1,nW,'Padding','compact','TileSpacing','compact');

        for w=1:nW
            nexttile; hold on
            % Assemble grouped data & positions
            allY = []; grp = []; gnames = {};
            for j=1:J
                y = Delta{w,j};
                allY = [allY; y(:)]; %#ok<AGROW>
                grp  = [grp; j*ones(numel(y),1)]; %#ok<AGROW>
            end
            for j=1:J, gnames{j} = labels{j}; end %#ok<AGROW>

            boxplot(allY, grp, ...
                'Labels', gnames, ...
                'Whisker', 1.5, ...
                'Symbol','', ...
                'Colors', cmap, ...
                'PlotStyle','traditional', ...
                'LabelOrientation','inline');
            
            % --- style fix here ---
            set(findobj(gca,'Tag','Median'),'LineWidth',1.5);
            set(findobj(gca,'Tag','Box'),'LineWidth',1.2);
            set(findobj(gca,'Tag','Outliers'),'Marker','none');
            set(findobj(gca,'Tag','Whisker'),'LineStyle','-');
            ylim([-1 1]);
            line(xlim,[0 0],'Color','k','LineWidth',1.5,'LineStyle','-');

            % Aesthetics
            title(sprintf('win %d–%d', S.match{w}{1}.idxRef(1), S.match{w}{1}.idxRef(end)));
            ylabel(sprintf('%s', met),'Interpreter','none');
            grid on; box on; ylim([-1 1]); %ylim(auto_ylim(allY));
            set(gca,'FontSize',12,'FontName','Calibri','XTickLabelRotation',20);
        end

        title(t, sprintf('%s — %s', upper(which_space), strrep(met,'_','\_')));
        xlabel(t, 'Datasets (per window)');
    end
end

function rng = auto_ylim(y)
    y = y(~isnan(y));
    if isempty(y), rng = [-1 1]; return; end
    q = quantile(y,[0.02 0.98]);
    pad = 0.05*range(q);
    if pad==0, pad=0.1*max(1,abs(q(2))); end
    rng = [q(1)-pad, q(2)+pad];
end
