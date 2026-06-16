function [Xb_n, y_n] = pull_node_pack2_base(Xcells, ymat, n)
% Collect all edges for node n (row-wise), excluding self.
mask = true(size(ymat)); mask((1:size(ymat,1))+ (n-1)*size(ymat,1)) = false; %#ok>IDIV
y_n = ymat(n, :)'; y_n(n) = [];  % row n, drop self
Xb_n = zeros(numel(y_n), numel(Xcells));
for c=1:numel(Xcells)
    v = Xcells{c}(n, :)'; v(n) = [];
    Xb_n(:,c) = v;
end
end
