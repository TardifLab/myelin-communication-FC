function varargout = INhandler(args, nin, defaults)
%INHANDLER Minimal optional-input helper retained for legacy functions.
if nargin < 2 || isempty(nin), nin = numel(args); end
if nargin < 3, defaults = {}; end
n = max(numel(defaults), nin);
out = defaults;
if numel(out) < n, out(end+1:n) = {[]}; end
for i = 1:min(numel(args), n)
    if ~isempty(args{i}), out{i} = args{i}; end
end
varargout = out(1:nargout);
if nargout == 1 && numel(out) == 1
    varargout{1} = out{1};
end
end
