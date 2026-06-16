function FC_struct = make_fc_struct(Yz_m, ylbl)
% MAKE_FC_STRUCT  Convert cell arrays of FC to struct array expected by blockwise_rsn_predict_fc.
%
% Inputs
%   Yz_m : 1xD cell, each NxN FC matrix
%   ylbl : 1xD cellstr, names for each FC matrix
%
% Output
%   FC_struct : 1xD struct('name', <char>, 'FC', <NxN double>)

assert(iscell(Yz_m) && iscell(ylbl), 'Yz_m and ylbl must be cell arrays.');
assert(numel(Yz_m)==numel(ylbl), 'Yz_m and ylbl must have same length.');

D = numel(Yz_m);
FC_struct = repmat(struct('name','', 'FC',[]), 1, D);
for d = 1:D
    A = Yz_m{d};
    if ~ismatrix(A) || size(A,1)~=size(A,2)
        error('FC matrix %d is not square.', d);
    end
    % Force symmetry (average with transpose) and zero diagonal
    A = (A + A.')/2;
    A(1:size(A,1)+1:end) = 0;
    FC_struct(d).name = char(string(ylbl{d}));
    FC_struct(d).FC   = A;
end
end
