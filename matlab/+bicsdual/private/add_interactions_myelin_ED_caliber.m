function [INT_models, INT_labels] = add_interactions_myelin_ED_caliber(RD_models, RD_labels, Base_both_X, caliber_names, ED_name)
% ADD_INTERACTIONS_MYELIN_ED_CALIBER
% For each RD model (spec = [R myelin..., D myelin...]), build a Level-3
% full spec that includes ONLY:
%   1) Within-CM: caliber_R .* (R myelin_i), caliber_D .* (D myelin_j)
%   2) With geometry: ED .* (R myelin_i) and ED .* (D myelin_j)
% and NO cross-CM myelin×myelin terms.
%
% Inputs
%   RD_models    : 1xM cell, each {1xK cell} of matrices [R..., D...]
%   RD_labels    : 1xM cellstr, label per RD model
%   Base_both_X  : 1xP cell of base predictors for BOTH (expected order:
%                  {binary_R, caliber_R, ED, binary_D, caliber_D})
%   caliber_names: 1x2 cellstr for labeling (default {'caliber_R','caliber_D'})
%   ED_name      : char label for ED (default 'ED')
%
% Outputs
%   INT_models   : 1xM cell, each {1x(K + KR + KD + KE) cell} full spec:
%                  [main effects spec, within-CM myelin×caliber, myelin×ED]
%   INT_labels   : 1xM cellstr (propagates RD_labels; interaction detail is implicit)

if nargin < 4 || isempty(caliber_names), caliber_names = {'caliber_R','caliber_D'}; end
if nargin < 5 || isempty(ED_name),       ED_name      = 'ED';                      end

% Pull base matrices (assumed positions)
try
    caliber_R = Base_both_X{2};
    ED        = Base_both_X{3};
    caliber_D = Base_both_X{5};
catch
    error('Base_both_X must at least contain {.., caliber_R(2), ED(3), .., caliber_D(5)} in that order.');
end

INT_models = cell(size(RD_models));
INT_labels = RD_labels;  % keep the same high-level label per model

for k = 1:numel(RD_models)
    spec = RD_models{k};
    K    = numel(spec);
    if K == 0
        INT_models{k} = spec;
        continue;
    end

    % Split into R and D halves (R first, then D)
    h = floor(K/2); if h==0, h = K; end
    Rspec = spec(1:h);
    Dspec = spec(h+1:end);

    % Build interactions (NO cross-CM myelin×myelin)
    INT = {};

    % (1) Within-CM: caliber × myelin (R)
    for i = 1:numel(Rspec)
        INT{end+1} = caliber_R .* Rspec{i}; %#ok<AGROW>
    end

    % (1) Within-CM: caliber × myelin (D)
    for j = 1:numel(Dspec)
        INT{end+1} = caliber_D .* Dspec{j}; %#ok<AGROW>
    end

    % (2) Myelin × ED (both CMs)
    for i = 1:numel(Rspec)
        INT{end+1} = ED .* Rspec{i}; %#ok<AGROW>
    end
    for j = 1:numel(Dspec)
        INT{end+1} = ED .* Dspec{j}; %#ok<AGROW>
    end

    % Full Level-3 spec: main effects + interactions (no cross-CM myelin×myelin)
    INT_models{k} = [spec, INT];
end
end
