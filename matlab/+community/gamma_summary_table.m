function T = gamma_summary_table(R, cfg)
%GAMMA_SUMMARY_TABLE Convert fine-sweep results to a human-readable table.

n = numel(R.gamma);
index = (1:n).';
gamma = double(R.gamma(:));
if isfield(R,'n_communities') && numel(R.n_communities)==n
    n_communities = double(R.n_communities(:));
else
    n_communities = arrayfun(@(k) numel(unique(R.consensus_partitions(:,k))), (1:n)).';
end
if isfield(R,'mean_modularity') && numel(R.mean_modularity)==n
    mean_modularity = double(R.mean_modularity(:));
else
    mean_modularity = mean(double(R.modularity),2,'omitnan');
end
zrand_mean = double(R.zrand_mean(:));
zrand_variance = double(R.zrand_variance(:));
selection_score = community.compute_selection_score(R, cfg.score_alpha);
is_score_selected = false(n,1);
if any(isfinite(selection_score))
    [~,idx] = max(selection_score,[],'omitnan');
    is_score_selected(idx) = true;
end
is_paper_selected = false(n,1);
if isfield(R,'paper') && isfield(R.paper,'selected_index') && ~isempty(R.paper.selected_index)
    idx = double(R.paper.selected_index);
    if idx>=1 && idx<=n, is_paper_selected(idx)=true; end
end
T = table(index,gamma,n_communities,mean_modularity,zrand_mean, ...
    zrand_variance,selection_score,is_score_selected,is_paper_selected);
end
