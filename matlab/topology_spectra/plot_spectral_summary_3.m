function plot_spectral_summary_3(results_D1,results_D2,results_D3,cmap,labels,str,plotWhich)
% Summarize spectral similarity metrics for three networks
%
% Input:
%   results_D1..D3 : outputs of spectral_alignment_compare_adj_lap()
%   cmap           : 3x3 colormap (rows correspond to D1..D3)
%   labels         : 1x3 cell of data labels (e.g., {'MTsat','Delay','Caliber'})
%   str            : char vector for figure title suffix
%   plotWhich      : optional 1x3 logical or numeric mask for line plot inclusion
%                    (default = [1 1 1]); true -> include in line plot
%
% Displays:
%   Left panel  : Cumulative cosine similarity curves (Adj solid, Lap dashed)
%                 Only datasets flagged by plotWhich are drawn (to reduce overlap).
%   Right panel : Wasserstein distances, grouped by operator (Adjacency vs Laplacian),
%                 with bars colored by the dataset row in cmap.
%
% 2025 Mark C Nelson (MNI)
%--------------------------------------------------------------------------

    if nargin < 7 || isempty(plotWhich)
        plotWhich = [1 1 1];
    end
    plotWhich = logical(plotWhich(:).'); % ensure 1x3 logical
    fntsz=20; fntname='Calibri';

    % --- Extract metrics
    cos_adj{1} = results_D1.cos_sim_adj;  cos_lap{1} = results_D1.cos_sim_lap;
    cos_adj{2} = results_D2.cos_sim_adj;  cos_lap{2} = results_D2.cos_sim_lap;
    cos_adj{3} = results_D3.cos_sim_adj;  cos_lap{3} = results_D3.cos_sim_lap;

    W_adj = [results_D1.wasserstein_adj, results_D2.wasserstein_adj, results_D3.wasserstein_adj];
    W_lap = [results_D1.wasserstein_lap, results_D2.wasserstein_lap, results_D3.wasserstein_lap];

    % --- Basic checks
    assert(all(cellfun(@(x) isnumeric(x)&&isvector(x), cos_adj)),'cos_sim_adj must be numeric vectors.');
    assert(all(cellfun(@(x) isnumeric(x)&&isvector(x), cos_lap)),'cos_sim_lap must be numeric vectors.');
    kmaxs = cellfun(@numel, cos_adj);
    assert(all(kmaxs == kmaxs(1)), 'All cos_sim vectors must have the same length.');
    kmax = kmaxs(1);
    assert(size(cmap,1) >= 3 && size(cmap,2) == 3, 'cmap must be 3x3 (RGB).');
    assert(numel(labels) == 3, 'labels must be a 1x3 cell array.');

    % --- Figure
    myfig(['Spectral comparison summary: ' str], [-2559 433 1102 364]); % [300 300 900 500]

    %% --- Panel 1: Cumulative Cosine Similarity (selectively plotted) ---
    subplot(1,2,1); cla; hold on;
    for d = 1:3
        if ~plotWhich(d), continue; end
        plot(1:kmax, cos_adj{d}, '-',  'LineWidth', 3, 'Color', cmap(d,:));
        plot(1:kmax, cos_lap{d}, '--', 'LineWidth', 3, 'Color', cmap(d,:));
    end
    hold off;
    xlabel('Leading eigenvectors (k)'); ylabel('Cumulative cosine similarity'); title('Spectral alignment');
    ylim([0 1]); xlim([0 kmax]); grid on; box on; font(fntsz,fntname);

    % Build legend only for included datasets (Adj/Lap entries per dataset)
    legtxt = {};
    for d = 1:3
        if plotWhich(d)
            legtxt{end+1} = [labels{d} ' - Adj']; %#ok<AGROW>
            legtxt{end+1} = [labels{d} ' - Lap']; %#ok<AGROW>
        end
    end
    if ~isempty(legtxt)
        legend(legtxt, 'Location','northeast');
    end

    %% --- Panel 2: Wasserstein summary (group Adj together, Lap together) ---
    % Data arranged as: rows = groups (Adj, Lap), columns = datasets (D1..D3)
    subplot(1,2,2); cla;
    data = [W_adj; W_lap];                    % 2 x 3
    B = bar(data, 'BarWidth', 0.8);           % creates 3 bar objects (one per dataset)
    % Color each dataset by its row in cmap
    for d = 1:3
        B(d).FaceColor = 'flat';
        B(d).CData = repmat(cmap(d,:), 2, 1); % same color for both groups
    end
    % X-axis: groups by operator
    set(gca,'XTick',1:2,'XTickLabel',{'Adjacency','Laplacian'});
    ylabel('Wasserstein distance'); title('Spectral distance');
    legend(labels, 'Location','northwest'); % legend indicates dataset colors
    ylim([0, max(data(:))*1.2 + eps]); grid on; set(gca,'XGrid','off'); font(fntsz,fntname);

    % Annotate values above each bar
    % bar() with grouped style positions bars around x = 1 (Adj) and x = 2 (Lap)
    % We can get XEndPoints from each Bar object. Fallback: compute offsets.
    try
        for d = 1:3
            xtips = B(d).XEndPoints;
            ytips = B(d).YEndPoints;
            for j = 1:numel(xtips)
                text(xtips(j), ytips(j) + 0.02*range(ylim), sprintf('%.2f', data(j,d)), ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontSize', 16, 'Color','k');
            end
        end
    catch
        % Fallback annotation: approximate offsets if XEndPoints not available
        xgrp = [1 2];
        nSeries = 3;
        w = 0.8;                      % BarWidth
        off = linspace(-w/2, w/2, nSeries+2);
        off = off(2:end-1);           % 3 offsets
        for j = 1:2
            for d = 1:3
                x = xgrp(j) + off(d);
                y = data(j,d);
                text(x, y + 0.02*range(ylim), sprintf('%.2f', y), ...
                    'HorizontalAlignment','center','VerticalAlignment','bottom', 'FontSize', 16, 'Color','k');
            end
        end
    end

end
