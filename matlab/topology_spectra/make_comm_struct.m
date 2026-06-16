function Comm_struct = make_comm_struct(useCM, useCMLBL, dsLabels, CMLBL_simple)
% MAKE_COMM_STRUCT  Flatten nested cells of communication matrices into a struct array.
%
% Inputs
%   useCM        : 1xS cell, each is 1xM cell of NxN communication matrices
%   useCMLBL     : 1xS cell, each is 1xM cell of labels combining structure+model (preferred)
%   dsLabels     : 1xS cell of structural labels (fallback if useCMLBL empty)
%   CMLBL_simple : 1xM cell of model labels (fallback if useCMLBL empty)
%
% Output
%   Comm_struct : 1x(S*M) struct with fields:
%       .name        : char (e.g., 'SPE_myelin')
%       .C           : NxN double (symmetrized)
%       .family      : 'path' or 'diff' (best-effort)
%       .weight_type : 'caliber'|'myelin'|'rate' (best-effort)

assert(iscell(useCM) && ~isempty(useCM), 'useCM must be a non-empty cell array.');

S = numel(useCM);
Comm_struct = struct('name',{},'C',{},'family',{},'weight_type',{});
k = 0;

for s = 1:S
    mats = useCM{s};
    assert(iscell(mats), 'useCM{%d} must be a 1xM cell of matrices.', s);
    M = numel(mats);

    % Decide labels to use
    if exist('useCMLBL','var') && ~isempty(useCMLBL) && ~isempty(useCMLBL{s})
        labels = useCMLBL{s};
        assert(numel(labels)==M, 'useCMLBL{%d} must match number of models.', s);
    else
        % Fallback: combine dsLabels + CMLBL_simple
        assert(iscell(dsLabels) && iscell(CMLBL_simple), 'Provide dsLabels and CMLBL_simple if useCMLBL is empty.');
        assert(numel(dsLabels)>=s, 'dsLabels shorter than useCM.');
        assert(numel(CMLBL_simple)==M, 'CMLBL_simple must match number of models.');
        labels = cell(1,M);
        for m = 1:M
            labels{m} = sprintf('%s_%s', char(string(CMLBL_simple{m})), char(string(dsLabels{s})));
        end
    end

    % Build struct entries
    for m = 1:M
        C = mats{m};
        if ~ismatrix(C) || size(C,1)~=size(C,2)
            error('Comm matrix at {%d}{%d} is not square.', s, m);
        end
        C = (C + C.')/2;                    % symmetrize
        C(1:size(C,1)+1:end) = 0;           % zero diagonal

        lbl = char(string(labels{m}));

        k = k+1;
        Comm_struct(k).name = lbl;
        Comm_struct(k).C    = C;

        % Best-effort family & weight_type from label text (tweak as needed)
        low = lower(lbl);
        if contains(low,'spe') || contains(low,'ne') || contains(low,'short') || contains(low,'spath') || contains(low,'nav')
            Comm_struct(k).family = 'path';
        elseif contains(low,'comm') || contains(low,'diff') || contains(low,'walk') || contains(low,'heat') || contains(low,'rw')
            Comm_struct(k).family = 'diff';
        else
            Comm_struct(k).family = 'unknown';
        end

        if contains(low,'myelin') || contains(low,'mtsat') || contains(low,'r1') || contains(low,'g-ratio') || contains(low,'gratio')
            Comm_struct(k).weight_type = 'myelin';
        elseif contains(low,'caliber') || contains(low,'diam') || contains(low,'length')
            Comm_struct(k).weight_type = 'caliber';
        elseif contains(low,'rate') || contains(low,'speed') || contains(low,'invdelay') || contains(low,'delay')
            % If you treat "delay" as a length-like cost for path models and 1/delay ("rate") as weight-like for diff models,
            % label them accordingly in useCMLBL for cleaner inference.
            if contains(low,'rate') || contains(low,'speed') || contains(low,'invdelay')
                Comm_struct(k).weight_type = 'rate';
            else
                Comm_struct(k).weight_type = 'caliber'; % neutral fallback; change if you prefer a distinct 'delay'
            end
        else
            Comm_struct(k).weight_type = 'unknown';
        end
    end
end
end
