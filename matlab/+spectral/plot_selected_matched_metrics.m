function figs = plot_selected_matched_metrics(T, config)
%PLOT_SELECTED_MATCHED_METRICS Plot selected matched-eigenvector metrics.
%
% The selected rows can be displayed in three ways:
%   'windowed' : supplemental-style grouped bars. Each tile is one metric,
%                x-axis groups are spectral windows, and bars are datasets.
%   'pooled'   : boxcharts pooling selected eigenvectors across windows.
%   'both'     : create both figure types.
%
% Relevant configuration fields:
%   config.selected_metrics_plot_mode : 'windowed' (default), 'pooled', 'both'
%   config.selected_metrics_bar_mode  : 'norm' (default) or 'raw'
%   config.selected_metric_names
%   config.figure_visible
%
% Windowed values are means across selected reference eigenvectors falling
% in each spectral window. In 'norm' mode, each metric is divided by its
% maximum absolute window mean within an operator, matching the legacy
% supplemental plotting convention.

if isempty(T)
    warning('No selected matched-eigenvector rows are available to plot.');
    figs = gobjects(0);
    return
end

plot_mode = local_config(config, 'selected_metrics_plot_mode', 'windowed');
bar_mode  = local_config(config, 'selected_metrics_bar_mode',  'norm');
figure_visible = local_config(config, 'figure_visible', 'on');
metrics = local_config(config, 'selected_metric_names', ...
    {'delta_SAr','delta_smooth','delta_PR','delta_modR2','delta_MoranI'});

plot_mode = lower(char(string(plot_mode)));
bar_mode  = lower(char(string(bar_mode)));

if ~ismember(plot_mode, {'windowed','pooled','both'})
    error('Unknown selected_metrics_plot_mode: %s', plot_mode);
end
if ~ismember(bar_mode, {'norm','raw'})
    error('Unknown selected_metrics_bar_mode: %s', bar_mode);
end

available_ops = cellstr(unique(string(T.operator), 'stable'));
operators = intersect({'adj','lap'}, available_ops, 'stable');
figs = gobjects(0);

for oi = 1:numel(operators)
    op = operators{oi};
    S = T(strcmp(string(T.operator), op), :);

    if isempty(S)
        continue
    end

    if ismember(plot_mode, {'windowed','both'})
        figs(end+1,1) = local_windowed_plot( ... %#ok<AGROW>
            S, op, metrics, bar_mode, figure_visible);
    end

    if ismember(plot_mode, {'pooled','both'})
        figs(end+1,1) = local_pooled_plot( ... %#ok<AGROW>
            S, op, metrics, figure_visible);
    end
end
end

% -------------------------------------------------------------------------
function fig = local_windowed_plot(S, op, metrics, bar_mode, figure_visible)
% Supplemental-style summary: one tile per metric, x-axis = windows.

datasets = unique(string(S.target_dataset), 'stable');
window_ids = unique(double(S.window_id));
window_ids = sort(window_ids(:));

n_datasets = numel(datasets);
n_windows = numel(window_ids);
n_metrics = numel(metrics);
colors = lines(n_datasets);

window_labels = strings(n_windows,1);
for wi = 1:n_windows
    rows = double(S.window_id) == window_ids(wi);
    kmin = min(double(S.window_k_min(rows)), [], 'omitnan');
    kmax = max(double(S.window_k_max(rows)), [], 'omitnan');
    window_labels(wi) = sprintf('k %d-%d', kmin, kmax);
end

fig = figure('Color','w', 'Visible',figure_visible, ...
    'Name',sprintf('Selected matched-eigenvector metrics windowed %s', ...
    upper(op)), 'NumberTitle','off');
layout = tiledlayout(1, n_metrics, ...
    'Padding','compact', 'TileSpacing','compact');

for mi = 1:n_metrics
    ax = nexttile(layout);
    metric = metrics{mi};

    if ~ismember(metric, S.Properties.VariableNames)
        axis(ax, 'off');
        title(ax, local_pretty_metric(metric));
        text(ax, 0.5, 0.5, 'Metric unavailable', ...
            'HorizontalAlignment','center');
        continue
    end

    % Rows = windows; columns = target datasets.
    values = nan(n_windows, n_datasets);
    for wi = 1:n_windows
        for di = 1:n_datasets
            keep = double(S.window_id) == window_ids(wi) & ...
                strcmp(string(S.target_dataset), datasets(di));
            x = double(S.(metric)(keep));
            x = x(isfinite(x));
            if ~isempty(x)
                values(wi,di) = mean(x);
            end
        end
    end

    scale_text = '';
    if strcmp(bar_mode, 'norm')
        finite_values = abs(values(isfinite(values)));
        if ~isempty(finite_values)
            scale = max(finite_values);
            if scale > 0
                values = values ./ scale;
            end
        end
        scale_text = ' (relative)';
    end

    if all(~isfinite(values(:)))
        axis(ax, 'off');
        title(ax, local_pretty_metric(metric));
        text(ax, 0.5, 0.5, 'No finite values', ...
            'HorizontalAlignment','center');
        continue
    end

    hold(ax, 'on');
    bars = bar(ax, 1:n_windows, values, 'grouped', 'BarWidth',0.78);
    for di = 1:min(numel(bars), n_datasets)
        bars(di).FaceColor = colors(di,:);
    end
    yline(ax, 0, 'k-', 'LineWidth',0.75);
    grid(ax, 'on');
    ax.XGrid = 'off';
    box(ax, 'on');
    ax.XTick = 1:n_windows;
    ax.XTickLabel = cellstr(window_labels);
    ax.XTickLabelRotation = 30;
    ax.TickLength = [0 0];
    xlim(ax, [0.5 n_windows + 0.5]);

    if strcmp(bar_mode, 'norm')
        ylim(ax, [-1 1]);
    end

    title(ax, local_pretty_metric(metric), 'Interpreter','tex');
    if mi == 1
        ylabel(ax, ['Target - caliber' scale_text]);
    end
    if mi == n_metrics
        legend(ax, bars, cellstr(datasets), ...
            'Location','northoutside', 'Orientation','horizontal');
    end
    hold(ax, 'off');
end

title(layout, sprintf('%s: selected matched eigenvectors by spectral window', ...
    upper(op)));
end

% -------------------------------------------------------------------------
function fig = local_pooled_plot(S, op, metrics, figure_visible)
% Original compact view: selected eigenvectors pooled across windows.

datasets = unique(string(S.target_dataset), 'stable');
n_datasets = numel(datasets);
colors = lines(n_datasets);

fig = figure('Color','w', 'Visible',figure_visible, ...
    'Name',sprintf('Selected matched-eigenvector metrics pooled %s', ...
    upper(op)), 'NumberTitle','off');
layout = tiledlayout(1, numel(metrics), ...
    'Padding','compact', 'TileSpacing','compact');

for mi = 1:numel(metrics)
    ax = nexttile(layout);
    metric = metrics{mi};

    if ~ismember(metric, S.Properties.VariableNames)
        axis(ax, 'off');
        title(ax, local_pretty_metric(metric));
        text(ax, 0.5, 0.5, 'Metric unavailable', ...
            'HorizontalAlignment','center');
        continue
    end

    hold(ax, 'on');
    plotted = false;
    for di = 1:n_datasets
        x = double(S.(metric)(strcmp(string(S.target_dataset), datasets(di))));
        x = x(isfinite(x));
        if isempty(x)
            continue
        end
        h = boxchart(ax, di * ones(size(x)), x);
        h.BoxFaceColor = colors(di,:);
        plotted = true;
    end

    if ~plotted
        axis(ax, 'off');
        title(ax, local_pretty_metric(metric));
        text(ax, 0.5, 0.5, 'No finite values', ...
            'HorizontalAlignment','center');
        continue
    end

    yline(ax, 0, 'k-', 'LineWidth',0.75);
    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [0.5 n_datasets + 0.5]);
    ax.XTick = 1:n_datasets;
    ax.XTickLabel = cellstr(datasets);
    ax.XTickLabelRotation = 30;
    title(ax, local_pretty_metric(metric), 'Interpreter','tex');
    if mi == 1
        ylabel(ax, 'Target - caliber');
    end
    hold(ax, 'off');
end

title(layout, sprintf('%s: selected matched eigenvectors pooled across windows', ...
    upper(op)));
end

% -------------------------------------------------------------------------
function value = local_config(config, field, default_value)
if isfield(config, field) && ~isempty(config.(field))
    value = config.(field);
else
    value = default_value;
end
end

function label = local_pretty_metric(metric)
switch lower(metric)
    case 'delta_sar'
        label = '\Delta S-A correlation';
    case 'delta_smooth'
        label = '\Delta smoothness';
    case 'delta_pr'
        label = '\Delta participation ratio';
    case 'delta_modr2'
        label = '\Delta module R^2';
    case 'delta_morani'
        label = '\Delta Moran''s I';
    otherwise
        label = strrep(metric, '_', '\_');
end
end
