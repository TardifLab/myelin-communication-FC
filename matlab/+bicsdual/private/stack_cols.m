function X = stack_cols(Xcells)
% Column-stack NxN cells into a  (nEdge x k) design.
if isempty(Xcells), X = zeros(0,0); return; end
v = cellfun(@vectorize_upper, Xcells, 'UniformOutput', false);
X = cat(2, v{:});
end
