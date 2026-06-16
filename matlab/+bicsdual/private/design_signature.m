function S = design_signature(X, names, kpeek)
% Quick per-column stats & stable hash for human inspection.
if nargin<3, kpeek = min(3, size(X,2)); end
if nargin<2 || isempty(names), names = arrayfun(@(i) sprintf('col%02d',i),1:size(X,2),'uni',0); end

% ---- Norm size & class or names ----
P = size(X, 2);

% normalize 'names' to 1xP cellstr
if nargin < 2 || isempty(names)
    names = arrayfun(@(j) sprintf('x%02d', j), 1:P, 'uni', 0);
elseif isstring(names)
    names = cellstr(names(:).'); % string -> cellstr row
elseif ischar(names)
    % Try to parse patterns like: BASE{a + b + c}
    s = strtrim(names);
    base = '';
    inner = '';
    m = regexp(s, '^(?<base>[^{}]*)\{(?<inner>[^}]*)\}\s*$', 'names');
    if ~isempty(m)
        base  = strtrim(m.base);
        inner = strtrim(m.inner);
        if ~isempty(inner)
            toks = regexp(inner, '\s*\+\s*', 'split');
            toks = cellfun(@strtrim, toks, 'uni', 0);
        else
            toks = {};
        end
        % If we parsed tokens and they match P, use them.
        if ~isempty(toks) && numel(toks) == P
            names = toks;
        elseif ~isempty(toks) && numel(toks) ~= P
            % Length mismatch: pad/trim to P
            if numel(toks) < P
                toks = [toks, arrayfun(@(j) sprintf('%s_tok%02d', base, j), 1:(P-numel(toks)), 'uni', 0)];
            else
                toks = toks(1:P);
            end
            names = toks;
        else
            % No tokens found: replicate base or generate defaults
            label = ~isempty(base) * base + isempty(base) * "x";
            names = arrayfun(@(j) sprintf('%s%02d', label, j), 1:P, 'uni', 0);
        end
    else
        % Plain single label: replicate or enumerate
        base = s;
        names = arrayfun(@(j) sprintf('%s%02d', base, j), 1:P, 'uni', 0);
    end
elseif iscell(names)
    % ensure row vector
    if iscolumn(names), names = names.'; end
    % If single cell with compound string, try parse like above
    if numel(names)==1 && (ischar(names{1}) || (isstring(names{1}) && isscalar(names{1})))
        tmp = char(string(names{1}));
        % recurse via char pathway
        s = strtrim(tmp);
        base = '';
        inner = '';
        m = regexp(s, '^(?<base>[^{}]*)\{(?<inner>[^}]*)\}\s*$', 'names');
        if ~isempty(m)
            base  = strtrim(m.base);
            inner = strtrim(m.inner);
            if ~isempty(inner)
                toks = regexp(inner, '\s*\+\s*', 'split');
                toks = cellfun(@strtrim, toks, 'uni', 0);
                if numel(toks) ~= P
                    if numel(toks) < P
                        toks = [toks, arrayfun(@(j) sprintf('%s_tok%02d', base, j), 1:(P-numel(toks)), 'uni', 0)];
                    else
                        toks = toks(1:P);
                    end
                end
                names = toks;
            else
                names = arrayfun(@(j) sprintf('%s%02d', base, j), 1:P, 'uni', 0);
            end
        else
            names = arrayfun(@(j) sprintf('%s%02d', s, j), 1:P, 'uni', 0);
        end
    end
else
    % last resort: try to stringify and enumerate
    names = arrayfun(@(j) sprintf('x%02d', j), 1:P, 'uni', 0);
end

% final guard: exact length
if numel(names) ~= P
    % pad or trim
    if numel(names) < P
        names = [names, arrayfun(@(j) sprintf('x%02d', j), 1:(P-numel(names)), 'uni', 0)];
    else
        names = names(1:P);
    end
end
% ---- end normalization block ----

S = struct('names',{names}, 'n',size(X,1), 'p',size(X,2));
S.col = struct('name',{},'min',{},'max',{},'mean',{},'std',{},'hash',{});
for j=1:size(X,2)
    x = X(:,j);
    S.col(j).name = names{j};
    S.col(j).min  = min(x);
    S.col(j).max  = max(x);
    S.col(j).mean = mean(x);
    S.col(j).std  = std(x);
    S.col(j).hash = md5_vec(x);  % below
end
% peek a few columns
S.peek_cols = {names{1:kpeek}};
end

function h = md5_vec(x)
% stable md5 on double vector
x = double(x(:));
b = typecast(x,'uint8');  % no I/O, purely deterministic in MATLAB session
% h = string(DataHash(b,'md5')); % requires FileExchange DataHash; fallback:
% If no DataHash, replace with: 
h = string(java.security.MessageDigest.getInstance('MD5').digest(uint8(b)));
end
