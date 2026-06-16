function paths = write_outputs(R, selected, cfg)
%WRITE_OUTPUTS Save summary tables, selected partition, and optional MAT/figures.

network_dir = fullfile(cfg.out_dir, lower(community.canonical_name(R.network_name)));
if ~exist(network_dir,'dir'), mkdir(network_dir); end
paths = struct();

if cfg.export_csv
    T = community.gamma_summary_table(R,cfg);
    paths.gamma_summary = fullfile(network_dir,'gamma_summary.csv');
    writetable(T,paths.gamma_summary);

    if isfield(selected,'partition') && ~isempty(selected.partition)
        node_id = (1:numel(selected.partition)).';
        raw_community = double(selected.raw_partition(:));
        community_id = double(selected.partition(:));
        included = isfinite(community_id);
        S = table(node_id,raw_community,community_id,included, ...
            'VariableNames',{'node_id','raw_community','community','included'});
        paths.selected_partition = fullfile(network_dir,'selected_partition.csv');
        writetable(S,paths.selected_partition);
    end
end

if cfg.save_results_mat
    paths.results_mat = fullfile(network_dir,'community_results.mat');
    community_results = R; %#ok<NASGU>
    selected_partition = selected; %#ok<NASGU>
    community_config = cfg; %#ok<NASGU>
    save(paths.results_mat,'community_results','selected_partition','community_config','-v7.3');
end
end
