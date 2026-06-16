function H = plot_stability(R, cfg, selected)
%PLOT_STABILITY Plot gamma-dependent community and stability diagnostics.

if nargin < 3, selected = struct(); end
T = community.gamma_summary_table(R,cfg);
H = figure('Color','w','Name',[R.network_name ' community stability']);
tl = tiledlayout(H,2,2,'TileSpacing','compact','Padding','compact');

nexttile(tl);
semilogx(T.gamma,T.n_communities,'LineWidth',1.5);
xlabel('Gamma'); ylabel('Number of communities'); title('Consensus partition size');

nexttile(tl);
semilogx(T.gamma,T.zrand_mean,'LineWidth',1.5);
xlabel('Gamma'); ylabel('Mean z-Rand'); title('Partition stability: mean');

nexttile(tl);
semilogx(T.gamma,T.zrand_variance,'LineWidth',1.5);
xlabel('Gamma'); ylabel('Variance of z-Rand'); title('Partition stability: variance');

nexttile(tl);
semilogx(T.gamma,T.selection_score,'LineWidth',1.5);
xlabel('Gamma'); ylabel('Selection score'); title(sprintf('Automated candidate (alpha = %.2f)',cfg.score_alpha));

if isfield(selected,'selected_gamma') && ~isempty(selected.selected_gamma)
    ax = findall(H,'Type','axes');
    for k = 1:numel(ax)
        xline(ax(k),selected.selected_gamma,'--','LineWidth',1.25);
    end
end
title(tl,sprintf('%s community-detection diagnostics',R.network_name), ...
    'Interpreter','none');
end
