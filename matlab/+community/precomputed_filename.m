function filename = precomputed_filename(network_name)
%PRECOMPUTED_FILENAME Return bundled file name for a structural network.

name = community.canonical_name(network_name);
filename = [lower(name) '_community_results.mat'];
end
