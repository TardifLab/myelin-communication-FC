function selected = select_partition(R, cfg)
%SELECT_PARTITION Select a paper, score-based, gamma-based, or no partition.

mode = lower(char(string(cfg.selection_mode)));
selected = struct('mode',mode,'selected_index',[], 'selected_gamma',[], ...
    'raw_partition',[], 'partition',[], 'n_communities',0, ...
    'excluded_nodes',[], 'selection_score',[]);

switch mode
    case 'none'
        return

    case 'paper'
        if ~isfield(R,'paper') || ~isfield(R.paper,'partition') || isempty(R.paper.partition)
            error('No bundled paper partition is available for %s.', R.network_name);
        end
        selected.selected_index = double(R.paper.selected_index);
        selected.selected_gamma = double(R.paper.selected_gamma);
        selected.raw_partition = double(R.paper.raw_partition(:));
        selected.partition = double(R.paper.partition(:));
        selected.n_communities = double(R.paper.n_communities);
        selected.excluded_nodes = double(R.paper.excluded_nodes(:));

    case 'score'
        score = community.compute_selection_score(R, cfg.score_alpha);
        if ~any(isfinite(score))
            error('No valid z-Rand values are available for score-based selection.');
        end
        [~,idx] = max(score,[],'omitnan');
        selected = from_index(R, idx, cfg, mode);
        selected.selection_score = score(idx);

    case 'gamma'
        if isempty(cfg.selected_gamma) || ~isscalar(cfg.selected_gamma)
            error('cfg.selected_gamma must be a scalar when selection_mode = ''gamma''.');
        end
        [~,idx] = min(abs(double(R.gamma(:)) - double(cfg.selected_gamma)));
        selected = from_index(R, idx, cfg, mode);

    otherwise
        error('Unknown cfg.selection_mode: %s', cfg.selection_mode);
end
end

function selected = from_index(R, idx, cfg, mode)
raw = double(R.consensus_partitions(:,idx));
part = community.relabel_partition(raw, cfg.min_community_size);
selected = struct();
selected.mode = mode;
selected.selected_index = idx;
selected.selected_gamma = double(R.gamma(idx));
selected.raw_partition = raw;
selected.partition = part;
selected.n_communities = numel(unique(part(isfinite(part))));
selected.excluded_nodes = find(~isfinite(part));
selected.selection_score = [];
end
