function [INT_models, INT_labels] = add_interactions_on_RDboth(RD_models, RD_labels, Base_both_X, caliber_names)
% Generate interactions:
%  - cross-CM: all pairwise myelin terms between the left (R) and right (D) halves of RD model
%  - within-CM: caliber_R × (R myelin terms), caliber_D × (D myelin terms)
% Assume RD_models entries are {R-terms..., D-terms...}
INT_models = cell(size(RD_models));
INT_labels = RD_labels; %#ok<NASGU>
% identify caliber matrices (we get them directly from Base_both_X positions)
caliber_R = Base_both_X{2};
caliber_D = Base_both_X{5};
for k=1:numel(RD_models)
    spec = RD_models{k};
    K = numel(spec);
    % Heuristic split: first half belongs to R, second half to D.
    h = floor(K/2); if h==0, h=K; end
    Rspec = spec(1:h); Dspec = spec(h+1:end);
    INT = {};
    % cross-CM interactions
    for i=1:numel(Rspec), for j=1:numel(Dspec)
        INT{end+1} = Rspec{i} .* Dspec{j}; %#ok<AGROW>
    end, end
    % within-CM caliber × myelin
    for i=1:numel(Rspec), INT{end+1} = caliber_R .* Rspec{i}; end %#ok<AGROW>
    for j=1:numel(Dspec), INT{end+1} = caliber_D .* Dspec{j}; end %#ok<AGROW>
    INT_models{k} = [spec, INT];
end
end
