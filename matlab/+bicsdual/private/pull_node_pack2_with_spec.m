function [Xb_n, Xf_n, y_n] = pull_node_pack2_with_spec(BASE_and_RED, SPEC, ymat, node_ix)
% Node-wise collector aligned to "mains + spec".
[Xb_n, y_n] = pull_node_pack2_base(BASE_and_RED, ymat, node_ix);
Xspec = pull_node_pack2_base(SPEC, ymat, node_ix);
Xf_n  = [Xb_n, Xspec];
end


% % % function [Xb,Xf,y_n,spec_cols] = pull_node_pack2_with_spec(BASE_and_RED, SPEC, ymat, n)
% % % cells = [BASE_and_RED(:); SPEC(:)]';
% % % K     = numel(cells);
% % % 
% % % Xi = cell(1,K);
% % % for k=1:K, A = cells{k}; Xi{k} = A(:,n); end
% % % Y = ymat(:,n);
% % % 
% % % M    = [Xi{:}];
% % % mask = ~(all(M==0,2) | isnan(Y));
% % % M    = M(mask,:);  Y = Y(mask);
% % % 
% % % keep = std(M,0,1,'omitnan') > 1e-12;
% % % M    = M(:, keep);
% % % 
% % % bw        = sum( keep(1:numel(BASE_and_RED)) );
% % % Xb        = M(:, 1:bw);
% % % Xf        = M;
% % % y_n       = Y;
% % % spec_cols = (bw+1):size(M,2);
% % % end
