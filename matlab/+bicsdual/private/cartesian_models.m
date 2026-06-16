function [RD_models, RD_labels] = cartesian_models(R_models, R_labels, D_models, D_labels)
RD_models = {};
RD_labels = {};
for r=1:numel(R_models)
    for d=1:numel(D_models)
        RD_models{end+1} = [R_models{r}, D_models{d}]; %#ok<AGROW>
        RD_labels{end+1} = sprintf('R{%s} + D{%s}', R_labels{r}, D_labels{d}); %#ok<AGROW>
    end
end
end