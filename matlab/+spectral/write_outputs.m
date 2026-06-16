function write_outputs(results,config)
%WRITE_OUTPUTS Write spectral CSV/MAT outputs and optional figures.

if ~exist(config.out_dir,'dir'), mkdir(config.out_dir); end

if config.write_csv
    if isfield(results,'alignment')
        outdir=fullfile(config.out_dir,'alignment');
        if ~exist(outdir,'dir'), mkdir(outdir); end
        targets=fieldnames(results.alignment);
        for j=1:numel(targets)
            target=targets{j};
            R=results.alignment.(target);
            for op={'adj','lap'}
                name=op{1};
                S=R.(name);
                K=numel(S.diagonal_similarity);
                T=table(repmat(string(target),K,1), ...
                    repmat(string(name),K,1),(1:K)', ...
                    S.diagonal_similarity,S.cumulative_similarity, ...
                    S.eigenvalues_ref,S.eigenvalues_target, ...
                    repmat(S.wasserstein,K,1), ...
                    'VariableNames',{'target_dataset','operator','k', ...
                    'diagonal_similarity','cumulative_similarity', ...
                    'caliber_eigenvalue','target_eigenvalue', ...
                    'wasserstein_distance'});
                writetable(T,fullfile(outdir,sprintf( ...
                    'alignment_caliber_%s_%s.csv',lower(target),name)));
            end
        end
    end

    if isfield(results,'reference_basis')
        B=results.reference_basis;
        Kadj=numel(B.adj.eigenvalues);
        Klap=numel(B.lap.eigenvalues);
        K=max(Kadj,Klap);
        adj=nan(K,1); lap=nan(K,1);
        adj(1:Kadj)=B.adj.eigenvalues;
        lap(1:Klap)=B.lap.eigenvalues;
        T=table((1:K)',adj,lap,'VariableNames', ...
            {'k','adjacency_eigenvalue','laplacian_eigenvalue'});
        writetable(T,fullfile(config.out_dir,'reference_basis_eigenvalues.csv'));
    end

    if isfield(results,'fingerprints')
        outdir=fullfile(config.out_dir,'fingerprints');
        if ~exist(outdir,'dir'), mkdir(outdir); end
        F=results.fingerprints;
        W=numel(F.wins);
        J=numel(F.DatasetLabels);
        M=numel(F.ModelLabels);
        rows=cell(W*J*M,12);
        r=0;
        for w=1:W
            for j=1:J
                for m=1:M
                    r=r+1;
                    rows(r,:)={w,min(F.wins{w}),max(F.wins{w}), ...
                        numel(F.wins{w}),F.DatasetLabels{j}, ...
                        F.ModelLabels{m},F.OpForModel{m}, ...
                        F.phi(w,j,m),F.phi_raw(w,j,m),F.dphi(w,j,m), ...
                        F.coverage(1,j,m),config.reference_dataset};
                end
            end
        end
        T=cell2table(rows,'VariableNames',{'window_id','k_min','k_max', ...
            'n_modes','dataset','communication_model','operator', ...
            'fraction','raw_fraction','delta_vs_reference','coverage', ...
            'reference_dataset'});
        writetable(T,fullfile(outdir,'spectral_fingerprints.csv'));
        Tw=table((1:W)',cellfun(@min,F.wins)',cellfun(@max,F.wins)', ...
            cellfun(@numel,F.wins)', ...
            'VariableNames',{'window_id','k_min','k_max','n_modes'});
        writetable(Tw,fullfile(outdir,'spectral_windows.csv'));
    end

    if isfield(results,'matched_table')
        outdir=fullfile(config.out_dir,'matched_eigenvectors');
        if ~exist(outdir,'dir'), mkdir(outdir); end
        writetable(results.matched_table, ...
            fullfile(outdir,'matched_eigenvector_metrics_all.csv'));
        writetable(results.matched_table(results.matched_table.is_selected,:), ...
            fullfile(outdir,'matched_eigenvector_metrics_selected.csv'));

        idx=results.selected_indices;
        Tadj=table(repmat("adj",numel(idx.adj),1),idx.adj(:), ...
            'VariableNames',{'operator','eigenvector_index'});
        Tlap=table(repmat("lap",numel(idx.lap),1),idx.lap(:), ...
            'VariableNames',{'operator','eigenvector_index'});
        writetable([Tadj;Tlap], ...
            fullfile(outdir,'selected_eigenvector_indices_used.csv'));
    end
end

if config.save_mat
    save_results=results;
    if isfield(save_results,'figure_handles')
        save_results=rmfield(save_results,'figure_handles');
    end
    if isfield(save_results,'matched') && isfield(save_results.matched,'figs')
        save_results.matched=rmfield(save_results.matched,'figs');
    end
    if isfield(save_results,'fingerprints') && isfield(save_results.fingerprints,'fig')
        save_results.fingerprints.fig=[];
    end
    save(fullfile(config.out_dir,'spectral_results.mat'), ...
        'save_results','config','-v7.3');
end

if config.save_figures && isfield(results,'figure_handles')
    outdir=fullfile(config.out_dir,'figures');
    if ~exist(outdir,'dir'), mkdir(outdir); end
    figs=results.figure_handles;
    for i=1:numel(figs)
        if ~isgraphics(figs(i),'figure'), continue; end
        name=get(figs(i),'Name');
        if isempty(name), name=sprintf('spectral_figure_%02d',i); end
        name=regexprep(lower(name),'[^a-z0-9]+','_');
        filename=fullfile(outdir,sprintf('%s_%02d.%s', ...
            name,i,config.figure_format));
        if exist('exportgraphics','file') == 2
            exportgraphics(figs(i),filename);
        else
            saveas(figs(i),filename);
        end
    end
end
end
