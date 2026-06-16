function figs = plot_pairwise_alignment(result, target_label, config)
%PLOT_PAIRWISE_ALIGNMENT Lightweight built-in plots for Module 1.

figs = gobjects(2,1);
figs(1) = figure('Color','w','Visible',config.figure_visible, ...
    'Name',['Spectral alignment: caliber vs ' target_label], ...
    'NumberTitle','off');
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

nexttile;
imagesc(result.adj.similarity_matrix);
axis image;
colorbar;
xlabel([target_label ' eigenvectors']);
ylabel('Caliber eigenvectors');
title('Normalized adjacency similarity');
clim([0 1]);

nexttile;
plot(result.adj.cumulative_similarity,'LineWidth',1.5);
grid on;
xlabel('Leading eigenvectors, k');
ylabel('Cumulative mean |cosine|');
ylim([0 1]);
title('Adjacency cumulative alignment');

nexttile;
imagesc(result.lap.similarity_matrix);
axis image;
colorbar;
xlabel([target_label ' eigenvectors']);
ylabel('Caliber eigenvectors');
title('Laplacian similarity');
clim([0 1]);

nexttile;
plot(result.lap.cumulative_similarity,'LineWidth',1.5);
grid on;
xlabel('Leading eigenvectors, k');
ylabel('Cumulative mean |cosine|');
ylim([0 1]);
title('Laplacian cumulative alignment');

figs(2) = figure('Color','w','Visible',config.figure_visible, ...
    'Name',['Spectral values: caliber vs ' target_label], ...
    'NumberTitle','off');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile;
hold on;
plot(result.adj.eigenvalues_ref,'LineWidth',1.5);
plot(result.adj.eigenvalues_target,'LineWidth',1.5);
grid on;
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
title(sprintf('Adjacency spectra; W=%.4g',result.adj.wasserstein));
legend({'caliber',target_label},'Interpreter','none','Location','best');

nexttile;
hold on;
plot(result.lap.eigenvalues_ref,'LineWidth',1.5);
plot(result.lap.eigenvalues_target,'LineWidth',1.5);
grid on;
xlabel('Eigenvalue index');
ylabel('Eigenvalue');
title(sprintf('Laplacian spectra; W=%.4g',result.lap.wasserstein));
legend({'caliber',target_label},'Interpreter','none','Location','best');
end
