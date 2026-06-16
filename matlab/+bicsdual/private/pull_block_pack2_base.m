function [Xb_blk, y_blk] = pull_block_pack2_base(Xcells, ymat, io, jo)
% Shared base: extract edges within (io,jo) block for all cells in Xcells.
idx = block_index(io,jo,size(ymat,1));
y_blk = ymat(idx);
Xb_blk = zeros(numel(y_blk), numel(Xcells));
for c=1:numel(Xcells), Xb_blk(:,c) = Xcells{c}(idx); end
end
