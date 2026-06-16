function [Xb_blk, Xf_blk, y_blk, spec_cols] = pull_block_pack2_with_spec(BASE_and_RED, SPEC, ymat, io, jo)
% Block-level collector. Returns spec_cols indexing into Xf_blk.

[Xb_blk, y_blk] = pull_block_pack2_base(BASE_and_RED, ymat, io, jo);
Xspec  = pull_block_pack2_base(SPEC, ymat, io, jo);   % same rows as block
Xf_blk = [Xb_blk, Xspec];
spec_cols = size(Xb_blk,2) + (1:size(Xspec,2));
end



% % % function [Xb,Xf,y_blk,spec_cols] = pull_block_pack2_with_spec(BASE_and_RED, SPEC, ymat, io, jo)
% % % % Build reduced (Base+RED) and full (Base+RED+SPEC) with a COMMON mask/drop.
% % % % Returns:
% % % %   Xb        : reduced matrix (Base+RED) after joint mask/drop
% % % %   Xf        : full matrix   (Base+RED+SPEC) after joint mask/drop
% % % %   y_blk     : vectorized block response after the same mask
% % % %   spec_cols : column indices in Xf that correspond to SPEC (post-drop)
% % % 
% % % cells = [BASE_and_RED(:); SPEC(:)]';
% % % K     = numel(cells);
% % % 
% % % Xi = cell(1,K);
% % % for k=1:K, A = cells{k}; Xi{k} = A(io,jo); end
% % % Y = ymat(io,jo);
% % % 
% % % % vectorize block
% % % if numel(io)==numel(jo) && all(io(:)==jo(:))
% % %     lt = tril(true(numel(io)));
% % %     for k=1:K, t = Xi{k}; Xi{k} = t(lt); end
% % %     Y = Y(lt);
% % % else
% % %     for k=1:K, Xi{k} = Xi{k}(:); end
% % %     Y = Y(:);
% % % end
% % % 
% % % M    = [Xi{:}];
% % % mask = ~(all(M==0,2) | isnan(Y));
% % % M    = M(mask,:);  Y = Y(mask);
% % % 
% % % % drop near-constant columns jointly
% % % keep = std(M,0,1,'omitnan') > 1e-12;
% % % M    = M(:, keep);
% % % 
% % % bw        = sum( keep(1:numel(BASE_and_RED)) );  % retained Base+RED width
% % % Xb        = M(:, 1:bw);
% % % Xf        = M;
% % % y_blk     = Y;
% % % spec_cols = (bw+1):size(M,2);
% % % end
