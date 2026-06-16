function [models, labels] = build_models_strict(Xlist, L, modes)
% Build exactly what user asked for in `modes`:
% modes can be char or cellstr; allowed: 'single','pairs','all'

    if ischar(modes) || isstring(modes), modes = {char(modes)}; end
    modes = lower(string(modes));
    want_single = any(modes=="single");
    want_pairs  = any(modes=="pairs");
    want_all    = any(modes=="all");

    K = numel(Xlist); models = {}; labels = {};

    % singles
    if want_single
        for i=1:K
            models{end+1} = {Xlist{i}}; %#ok<AGROW>
            labels{end+1} = L{i};       %#ok<AGROW>
        end
    end

    % pairs
    if want_pairs && K>=2
        for i=1:K-1
            for j=i+1:K
                models{end+1} = {Xlist{i}, Xlist{j}}; %#ok<AGROW>
                labels{end+1} = [L{i} '+' L{j}];      %#ok<AGROW>
            end
        end
    end

    % all
    if want_all && K>=1
        models{end+1} = Xlist;         %#ok<AGROW>
        labels{end+1} = strjoin(L,'+'); %#ok<AGROW>
    end
end
