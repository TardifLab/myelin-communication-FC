function [partition, meta] = load_paper_partition(network_name, precomputed_dir)
%LOAD_PAPER_PARTITION Load the exact selected partition used in the paper.

if nargin < 2 || isempty(precomputed_dir)
    rootdir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    precomputed_dir = fullfile(rootdir,'data','community_detection','precomputed');
end
R = community.load_precomputed(network_name,precomputed_dir);
partition = double(R.paper.partition(:));
meta = R.paper;
meta.network_name = R.network_name;
end
