function Y = mylog(X, offset)
%MYLOG Safe natural log transform used by legacy code.
if nargin < 2 || isempty(offset), offset = 0; end
Y = X;
mask = isfinite(X) & X > 0;
Y(mask) = log(X(mask) + offset);
Y(~mask) = 0;
end
