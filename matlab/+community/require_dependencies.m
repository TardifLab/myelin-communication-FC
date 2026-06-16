function require_dependencies()
%REQUIRE_DEPENDENCIES Check functions required to rerun community detection.

required = {'community_louvain','agreement','consensus_und','fcn_randz'};
missing = required(cellfun(@(f) exist(f,'file') ~= 2, required));
if ~isempty(missing)
    error(['Missing community-detection dependencies: %s. ', ...
        'Add the Brain Connectivity Toolbox and the function fcn_randz to the MATLAB path.'], ...
        strjoin(missing, ', '));
end
end
