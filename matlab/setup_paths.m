function setup_paths(rootdir)
%SETUP_PATHS Add this repository's MATLAB code to the path.
%
% Usage:
%   setup_paths
%   setup_paths('/path/to/repository')

if nargin < 1 || isempty(rootdir)
    rootdir = fileparts(fileparts(mfilename('fullpath')));
end
addpath(genpath(fullfile(rootdir,'matlab')));
end
