function results=spectral_alignment_compare_adj_lap(A1,A2,kmax,doSignAlign,str)
% Compare spectral alignment between two weighted brain networks
% for both adjacency and Laplacian representations.
%
% Inputs:
%   A1, A2        : weighted adjacency matrices (NxN)
%   kmax          : number of leading eigenvectors to analyze
%   doSignAlign   : logical flag (true/false) to align eigenvector signs
%   str           : label for plots
%
% Outputs:
%   results : structure containing
%       .cos_sim_adj, .cos_sim_lap        - cumulative cosine similarities
%       .eigvec_corr_adj, .eigvec_corr_lap- full similarity matrices
%       .eigvals_adj1, .eigvals_adj2      - adjacency eigenvalue spectra
%       .eigvals_lap1, .eigvals_lap2      - Laplacian eigenvalue spectra
%       .wasserstein_adj, .wasserstein_lap- Wasserstein distances
%
% Figures: 
%   (1) 4-panel alignment:
%       [1] Adjacency eigenvector alignment heatmap
%       [2] Adjacency cumulative similarity curve
%       [3] Laplacian eigenvector alignment heatmap
%       [4] Laplacian cumulative similarity curve
%   (2) Joint spectral plot (Adjacency + Laplacian spectra overlayed)
%
% 2025 Mark C Nelson (MNI)
%--------------------------------------------------------------------------

    %% --- 1. Adjacency-based spectral alignment ---
    fprintf('--- Adjacency alignment ---\n');
    [cos_sim_adj, eigvec_corr_adj, eigvals_adj1, eigvals_adj2] = ...
        compute_alignment(A1, A2, kmax, doSignAlign, 'adjacency');

    %% --- 2. Laplacian-based spectral alignment ---
    fprintf('--- Laplacian alignment ---\n');
    L1 = diag(sum(A1,2)) - A1;
    L2 = diag(sum(A2,2)) - A2;
    [cos_sim_lap, eigvec_corr_lap, eigvals_lap1, eigvals_lap2] = ...
        compute_alignment(L1, L2, kmax, doSignAlign, 'laplacian');

    %% --- 3. Compute Wasserstein distances (spectral distance)
    wasserstein_adj = wasserstein_1D(eigvals_adj1, eigvals_adj2);
    wasserstein_lap = wasserstein_1D(eigvals_lap1, eigvals_lap2);

    %% --- 4. Plot results ---
    myfig(['Spectral alignment (adjacency & Laplacian): ' str],[-2559 -85 1233 882]);

    % Adjacency heatmap
    subplot(2,2,1); imagesc(eigvec_corr_adj);
    axis square; colorbar; colormap(smartcmaps('inferno')); clim([0 1]);
    xlabel('A2 eigenvectors'); ylabel('A1 eigenvectors');
    title('Adjacency eigenvector similarity'); font(20,'Cambria');

    % Adjacency cumulative
    subplot(2,2,2); plot(1:kmax,cos_sim_adj,'k-','LineWidth',3);
    xlabel('k (leading eigenvectors)'); ylabel('Mean cosine similarity');
    ylim([0 1]); grid on; font(20,'Cambria');
    title('Adjacency cumulative alignment');

    % Laplacian heatmap
    subplot(2,2,3); imagesc(eigvec_corr_lap);
    axis square; colorbar; colormap(smartcmaps('inferno')); clim([0 1]);
    xlabel('L2 eigenvectors'); ylabel('L1 eigenvectors');
    title('Laplacian eigenvector similarity'); font(20,'Cambria');

    % Laplacian cumulative
    subplot(2,2,4); plot(1:kmax,cos_sim_lap,'k-','LineWidth',3);
    xlabel('k (leading eigenvectors)'); ylabel('Mean cosine similarity');
    ylim([0 1]); grid on; font(20,'Cambria');
    title('Laplacian cumulative alignment');

    % sgtitle('Spectral alignment: adjacency vs Laplacian');

  %% --- 5. Joint spectral plot ---
    myfig(['Joint spectral comparison: ' str],[300 300 800 400]);

    subplot(1,2,1); hold on;
    plot(eigvals_adj1,'k.-','LineWidth',2,'DisplayName','Adjac - A1');
    plot(eigvals_adj2,'r.-','LineWidth',2,'DisplayName','Adjac - A2');
    xlabel('Eigenvalue index'); ylabel('Eigenvalue'); font(20,'Cambria');
    title(sprintf('A-spectra (W = %.2f)', wasserstein_adj));
    legend('show','Location','best');
    grid on; box off;

    subplot(1,2,2); 
    semilogy(eigvals_lap1,'k.-','LineWidth',2,'DisplayName','Laplac - L1'); 
    hold on;
    semilogy(eigvals_lap2,'r.-','LineWidth',2,'DisplayName','Laplac - L2');
    xlabel('Eigenvalue index'); ylabel('Eigenvalue (log scale)'); font(20,'Cambria');
    title(sprintf('L-spectra (W = %.2f)', wasserstein_lap));
    legend('show','Location','best');
    grid on; box off;

    % sgtitle('Joint spectral comparison');

    %% --- 5. Save outputs in struct ---
    results = struct( ...
        'cos_sim_adj', cos_sim_adj, ...
        'cos_sim_lap', cos_sim_lap, ...
        'eigvec_corr_adj', eigvec_corr_adj, ...
        'eigvec_corr_lap', eigvec_corr_lap, ...
        'eigvals_adj1', eigvals_adj1, ...
        'eigvals_adj2', eigvals_adj2, ...
        'eigvals_lap1', eigvals_lap1, ...
        'eigvals_lap2', eigvals_lap2, ...
        'wasserstein_adj', wasserstein_adj, ...
        'wasserstein_lap', wasserstein_lap ...
    );
end


%% ---------- Helper: Eigenvector alignment ----------
function [cos_sim_all,eigvec_corr_mat,eigvals1,eigvals2]=compute_alignment(M1,M2,kmax,doSignAlign,mode)
    % Ensure symmetry
    M1 = (M1 + M1') / 2;
    M2 = (M2 + M2') / 2;

    % Eigen decomposition
    switch lower(mode)
        case 'adjacency'
            [V1, D1] = eigs(M1, kmax, 'largestabs');
            [V2, D2] = eigs(M2, kmax, 'largestabs');

            eigvals1 = diag(D1);
            eigvals2 = diag(D2);

            % sort eigenpairs by |eigenvalue| to match 'largestabs'
            [~, idx1] = sort(abs(eigvals1), 'descend');
            [~, idx2] = sort(abs(eigvals2), 'descend');

        case 'laplacian'
            [V1, D1] = eigs(M1, kmax+1, 'smallestabs');
            [V2, D2] = eigs(M2, kmax+1, 'smallestabs');
            
            eigvals1 = diag(D1);
            eigvals2 = diag(D2);
    
            % sort eigenpairs from low to high frequency
            [eigvals1, idx1] = sort(eigvals1, 'ascend');
            [eigvals2, idx2] = sort(eigvals2, 'ascend');

        otherwise
            error('Unknown mode: choose ''adjacency'' or ''laplacian''');
    end

    % reorder eigenvectors
    V1 = V1(:, idx1);
    V2 = V2(:, idx2);
    
    % for laplacian, drop trivial eigenmode (assumes only 1 i.e., graph is connected)
    if strcmpi(mode, 'laplacian')
        eigvals1 = eigvals1(2:kmax+1);
        eigvals2 = eigvals2(2:kmax+1);
        V1 = V1(:, 2:kmax+1);
        V2 = V2(:, 2:kmax+1);
    end

    % for adjacency, also reorder eigenvalues after sorting
    if strcmpi(mode, 'adjacency')
        eigvals1 = eigvals1(idx1);
        eigvals2 = eigvals2(idx2);
    end

    % Normalize eigenvectors
    V1 = bsxfun(@rdivide, V1, sqrt(sum(V1.^2, 1)));
    V2 = bsxfun(@rdivide, V2, sqrt(sum(V2.^2, 1)));

    % Optional sign alignment
    if doSignAlign
        for i = 1:kmax
            % if corr(V1(:,i), V2(:,i)) < 0
            if dot(V1(:,i), V2(:,i)) < 0
                V2(:,i) = -V2(:,i);
            end
        end
    end

    % Cosine similarity matrix
    eigvec_corr_mat = abs(V1' * V2);
    eigvec_corr_mat = min(eigvec_corr_mat, 1);

    % Diagonal and cumulative means
    diag_corr = diag(eigvec_corr_mat);
    cos_sim_all = arrayfun(@(kk) mean(diag_corr(1:kk)), 1:kmax);
end


%% ---------- Helper: Wasserstein distance (1D spectra) ----------
function W = wasserstein_1D(lambda1, lambda2)

  % If weights equal (normalized, spectra signed)
    n = min(numel(lambda1), numel(lambda2));
    l1 = sort(lambda1(1:n));
    l2 = sort(lambda2(1:n));
    W  = mean(abs(l1 - l2));

  % % % % CDF of sums trick (not appropriate for signed spectra)
  % % % % Normalize length
  % % %   n = min(length(lambda1), length(lambda2));
  % % %   lambda1 = sort(lambda1(1:n));
  % % %   lambda2 = sort(lambda2(1:n));
  % % % 
  % % % % Compute empirical CDFs and their L1 difference
  % % %   F1 = cumsum(lambda1) / sum(lambda1);
  % % %   F2 = cumsum(lambda2) / sum(lambda2);
  % % %   W = mean(abs(F1 - F2));
end