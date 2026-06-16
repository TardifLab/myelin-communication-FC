function idx = block_index(io, jo, n)
% Linear indices for edges within RSN block (io,jo), i<j only.
mask = false(n);
mask(io, jo) = true; mask(jo, io) = true;
mask(1:n+1:end) = false;      % drop diagonal
mask = triu(mask,1);          % i<j
idx = find(mask);
end
