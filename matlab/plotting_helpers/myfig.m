function h = myfig(name, position)
%MYFIG Lightweight compatibility wrapper around MATLAB figure.
%
% This preserves compatibility with the original topology_spectra helpers
% while using only built-in graphics.

if nargin < 1 || isempty(name), name = ''; end
if nargin < 2 || isempty(position), position = [100 100 900 600]; end
h = figure('Color','w','Name',char(string(name)), ...
    'NumberTitle','off','Position',position);
end
