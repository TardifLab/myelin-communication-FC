function X = load_matrix_var(filename, varname)
%LOAD_MATRIX_VAR Load a square numeric matrix from .mat/.csv/.txt.
%
% For .mat files, this function first tries `varname`, then `Dtg`, then the
% first square numeric variable in the file.

if nargin < 2, varname = ''; end
[~,~,ext] = fileparts(filename);
if ~exist(filename,'file')
    error('Input file not found: %s', filename);
end
switch lower(ext)
    case '.mat'
        S = load(filename);
        if ~isempty(varname) && isfield(S,varname)
            X = S.(varname);
        elseif isfield(S,'Dtg')
            X = S.Dtg;
        else
            names = fieldnames(S);
            X = [];
            for i = 1:numel(names)
                candidate = S.(names{i});
                if isnumeric(candidate) && ismatrix(candidate) && size(candidate,1)==size(candidate,2)
                    X = candidate;
                    break
                end
            end
            if isempty(X)
                error('No square numeric matrix found in %s', filename);
            end
        end
    case {'.csv','.txt'}
        X = readmatrix(filename);
    otherwise
        error('Unsupported matrix file extension: %s', ext);
end
X = double(X);
if size(X,1) ~= size(X,2)
    error('Matrix must be square: %s', filename);
end
X(~isfinite(X)) = 0;
X(1:size(X,1)+1:end) = 0;
end
