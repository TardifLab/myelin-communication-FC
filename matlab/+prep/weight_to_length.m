function L = weight_to_length(W, method)
%WEIGHT_TO_LENGTH Convert weights to path lengths.
if nargin < 2 || isempty(method), method = 'log'; end
W = double(W);
W(~isfinite(W)) = 0;
switch lower(method)
    case 'log'
        if any(W(:) < 0 | W(:) > 1)
            error('For -log(w), weights must be normalized into [0,1].');
        end
        L = -log(W);
        L(isinf(L)) = 0;
    case 'inv'
        L = 1 ./ W;
        L(isinf(L)) = 0;
    otherwise
        error('Unknown length transform: %s', method);
end
L(~isfinite(L)) = 0;
mx = max(L(:));
if mx > 0, L = L ./ mx; end
L = (L + L.') ./ 2;
L(1:size(L,1)+1:end) = 0;
end
