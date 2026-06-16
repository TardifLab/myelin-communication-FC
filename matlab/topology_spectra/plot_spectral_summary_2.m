function plot_spectral_summary_2(results_D1,results_D2,cmap,labels,str)
% Summarize spectral similarity metrics for two networks
%
% Input:
%   results_D1, results_D2  : outputs of spectral_alignment_compare_adj_lap()
%   cmap                    : 2x3 colormap
%   labels                  : 1x2 cell of data labels
%   str                     : char vector for plot labels
%   
%
% Displays:
%   Generates one compact figure summarizing:
%     - Cumulative cosine similarity curves for adjacency & Laplacian
%     - Wasserstein distances as annotated bars
%     - Combined radar-style or overlay view to visualize "routing vs diffusion" overlap
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------


    % Extract metrics
    cos_adj_D1 = results_D1.cos_sim_adj;
    cos_lap_D1 = results_D1.cos_sim_lap;
    cos_adj_D2 = results_D2.cos_sim_adj;
    cos_lap_D2 = results_D2.cos_sim_lap;

    W_adj_D1 = results_D1.wasserstein_adj;
    W_lap_D1 = results_D1.wasserstein_lap;
    W_adj_D2 = results_D2.wasserstein_adj;
    W_lap_D2 = results_D2.wasserstein_lap;

    kmax = numel(cos_adj_D1);

    myfig(['Spectral comparison summary: ' str],[-2559 282 683 515]); % [300 300 900 500]

    %% --- Panel 1: Cumulative Cosine Similarity ---
    subplot(1,2,1)
    hold on;
    plot(1:kmax, cos_adj_D1,'-', 'LineWidth',3,'Color',cmap(1,:));
    plot(1:kmax, cos_lap_D1,'--','LineWidth',3,'Color',cmap(1,:));
    plot(1:kmax, cos_adj_D2,'-', 'LineWidth',3,'Color',cmap(2,:));
    plot(1:kmax, cos_lap_D2,'--','LineWidth',3,'Color',cmap(2,:));
    hold off;

    xlabel('Leading eigenvectors (k)'); ylabel('Mean cosine similarity'); title('Spectral alignment');
    ylim([0 1]); xlim([0 kmax]); grid on; box on; font(20,'Cambria');
    legend({[labels{1} '-Adj'],[labels{1} '-Lap'],[labels{2} '-Adj'],[labels{2} '-Lap']}, 'Location','northeast');
    

    % % % text(kmax*0.6, 0.9, 'Adj = solid', 'Color', [0.3 0.3 0.3]);
    % % % text(kmax*0.6, 0.85, 'Lap = dashed', 'Color', [0.3 0.3 0.3]);

    %% --- Panel 2: Wasserstein summary (bar chart) ---
    subplot(1,2,2)
    data = [W_adj_D1, W_lap_D1; W_adj_D2, W_lap_D2];
    B=bar(data,'FaceColor','flat','BarWidth',0.8);
    % B(1).CData=cmap; B(2).CData=cmap; B(2).FaceAlpha=0.3;
    B(1).CData=[.1 .1 .1; .1 .1 .1]; B(2).CData=[.7 .7 .7; .7 .7 .7];
    set(gca,'XTickLabel',labels);
    ylabel('Wasserstein distance'); title('Spectral distance');
    legend({'Adjacency','Laplacian'}, 'Location','northwest');
    ylim([0 0.6]); grid on; font(20,'Cambria');
    

    % Annotate values on bars
    for i = 1:2
        text(i-0.3,data(i,1)+0.02,sprintf('%.2f',data(i,1)),'FontSize',18,'Color','k');
        text(i+0.02,data(i,2)+0.02,sprintf('%.2f',data(i,2)),'FontSize',18,'Color','k');
    end

    % sgtitle('Spectral comparison summary: Routing vs Diffusion','FontSize',14,'FontWeight','bold');

end
