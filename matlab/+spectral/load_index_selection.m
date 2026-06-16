function indices = load_index_selection(config,max_k)
%LOAD_INDEX_SELECTION Load paper/custom/all matched eigenvector indices.

mode = lower(string(config.matched_indices_mode));
switch mode
    case "paper"
        if ~exist(config.matched_indices_csv,'file')
            error('Matched-eigenvector index CSV not found: %s',config.matched_indices_csv);
        end
        T = readtable(config.matched_indices_csv);
        required = {'eigenvector_index','use_adjacency','use_laplacian'};
        if ~all(ismember(required,T.Properties.VariableNames))
            error('Matched-index CSV is missing required columns.');
        end
        k = double(T.eigenvector_index);
        indices.adj = unique(k(logical(T.use_adjacency)),'stable');
        indices.lap = unique(k(logical(T.use_laplacian)),'stable');
    case "custom"
        indices = config.custom_indices;
    case "all"
        indices.adj = 1:max_k;
        indices.lap = 1:max_k;
    otherwise
        error('Unknown matched_indices_mode: %s',mode);
end
indices.adj = double(indices.adj(:)');
indices.lap = double(indices.lap(:)');
indices.adj = indices.adj(indices.adj>=1 & indices.adj<=max_k);
indices.lap = indices.lap(indices.lap>=1 & indices.lap<=max_k);
end
