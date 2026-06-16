function INT_models = apply_interactions_flagged(models_main, side, interaction, Base_both_X)
% APPEND interactions to a single "all-myelin" spec for one CM side or both.
% models_main : 1xM cell, each {1xK cell} (we assume M==1 when modes='all')
% side        : 'R' | 'D' | 'RD'   (R-only, D-only, or both-halves present)
% interaction : 'none' | 'caliber' | 'ed' | 'both'
% Base_both_X : EITHER:
%               - 1x5 both-sides pack: {bin_R, cal_R, ED, bin_D, cal_D}
%               - 1x3 side-specific  : {bin_S, cal_S, ED}  (S determined by `side`)
%
% Notes:
% - Robust to 1x3 or 1x5 inputs.
% - If side='RD', a 1x5 pack is required (throws a clear error otherwise).

if strcmpi(interaction,'none')
    INT_models = models_main; return;
end

% ----------------------------
% Normalize bases by shape
% ----------------------------
n = numel(Base_both_X);
switch n
    case 5
        % Expected both-sides order
        bin_R = Base_both_X{1};
        cal_R = Base_both_X{2};
        ED    = Base_both_X{3};
        bin_D = Base_both_X{4};
        cal_D = Base_both_X{5};

    case 3
        % Side-specific trio: {bin_S, cal_S, ED}
        % Map to R or D depending on `side`.
        ED = Base_both_X{3};
        switch upper(side)
            case 'R'
                bin_R = Base_both_X{1};
                cal_R = Base_both_X{2};
                bin_D = [];
                cal_D = [];
            case 'D'
                bin_R = [];
                cal_R = [];
                bin_D = Base_both_X{1};
                cal_D = Base_both_X{2};
            case 'RD'
                error('apply_interactions_flagged:BothSidesRequired', ...
                    ['side="RD" requires Base_both_X as 1x5: {bin_R, cal_R, ED, bin_D, cal_D}. ', ...
                     'Received a 1x3 side-specific pack instead.']);
            otherwise
                error('Unknown side: %s', side);
        end

    otherwise
        error('apply_interactions_flagged:BadShape', ...
              'Base_both_X must have 3 or 5 elements. Got %d.', n);
end

if isempty(ED)
    error('apply_interactions_flagged:MissingED', 'ED must be provided (element 3).');
end

% ----------------------------
% Build interactions
% ----------------------------
INT_models = cell(size(models_main));
for k = 1:numel(models_main)
    spec = models_main{k};              % {R..., D...} or just one side
    K    = numel(spec);
    if K == 0
        INT_models{k} = spec;
        continue;
    end

    % Split spec by requested side
    switch upper(side)
        case 'R'
            Rspec = spec;   Dspec = {};
        case 'D'
            Rspec = {};     Dspec = spec;
        case 'RD'
            h     = floor(K/2);
            if h == 0, h = K; end
            Rspec = spec(1:h);
            Dspec = spec(h+1:end);
        otherwise
            error('Unknown side: %s', side);
    end

    INT = {};

    % caliber interactions
    if any(strcmpi(interaction, {'caliber','both'}))
        if ~isempty(cal_R) && ~isempty(Rspec)
            for i = 1:numel(Rspec), INT{end+1} = cal_R .* Rspec{i}; end %#ok<AGROW>
        end
        if ~isempty(cal_D) && ~isempty(Dspec)
            for j = 1:numel(Dspec), INT{end+1} = cal_D .* Dspec{j}; end %#ok<AGROW>
        end
    end

    % ED interactions (shared ED for both sides)
    if any(strcmpi(interaction, {'ed','both'}))
        if ~isempty(Rspec)
            for i = 1:numel(Rspec), INT{end+1} = ED .* Rspec{i}; end %#ok<AGROW>
        end
        if ~isempty(Dspec)
            for j = 1:numel(Dspec), INT{end+1} = ED .* Dspec{j}; end %#ok<AGROW>
        end
    end

    % append interactions
    INT_models{k} = [spec, INT];
end
end