function [models, labels] = build_models(Xlist, L, doPairs, doAll)
if nargin<3, doPairs=true; end
if nargin<4, doAll=true; end
K = numel(Xlist); models = {}; labels = {};
% singles
for i=1:K, models{end+1} = {Xlist{i}}; labels{end+1} = L{i}; end %#ok<AGROW>
% pairs
if doPairs
    for i=1:K-1, for j=i+1:K
        models{end+1} = {Xlist{i}, Xlist{j}}; labels{end+1} = [L{i} '+' L{j}];
    end, end
end
% all
if doAll && K>2
    models{end+1} = Xlist; labels{end+1} = strjoin(L,'+');
end
end