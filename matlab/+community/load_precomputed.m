function R = load_precomputed(network_name, precomputed_dir)
%LOAD_PRECOMPUTED Load and validate a standardized gamma-sweep result.

name = community.canonical_name(network_name);
filename = fullfile(precomputed_dir, community.precomputed_filename(name));
if ~exist(filename, 'file')
    error('Precomputed community file not found: %s', filename);
end

S = load(filename);
if ~isfield(S, 'community_results')
    error('Expected variable community_results in %s.', filename);
end
R = S.community_results;

required = {'network_name','gamma','consensus_partitions','zrand_mean', ...
    'zrand_variance','modularity','input_matrix','transformed_matrix','paper'};
for k = 1:numel(required)
    if ~isfield(R, required{k})
        error('Missing field community_results.%s in %s.', required{k}, filename);
    end
end

R.network_name = name;
R.source_file = filename;
R.mode = 'precomputed';
end
