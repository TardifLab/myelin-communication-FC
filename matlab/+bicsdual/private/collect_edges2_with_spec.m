function [Xb, Xf, y] = collect_edges2_with_spec(BASE_and_RED, SPEC, ymat)
% Collect all upper (or full) edges into vector form for GLOBAL ΔR².
% Left block (base side) = [BASE RED], right block (full) = [BASE RED SPEC].
% Uses an intercept INSIDE safe_compare3 (not here).

% Vectorize y
y = vectorize_upper(ymat);

% Concatenate columns in the given order
Xb = stack_cols(BASE_and_RED);
Xf = [Xb, stack_cols(SPEC)];
end



% % % function [Xb,Xf,y_vec,spec_cols] = collect_edges2_with_spec(BASE_and_RED, SPEC, ymat)
% % % % Vectorize full matrix, build joint mask/drop over [BASE RED SPEC], then split.
% % % % BASE_and_RED : 1×B cell (BASE then RED mains)
% % % % SPEC         : 1×S cell (e.g., interactions only)
% % % % ymat         : NxN FC matrix
% % % 
% % % N = size(ymat,1);
% % % cells = [BASE_and_RED(:); SPEC(:)]';
% % % K     = numel(cells);
% % % 
% % % % vectorize all edges (use lower-tri incl diag to match your pipeline if preferred)
% % % lt = tril(true(N));
% % % Xi = cell(1,K);
% % % for k=1:K
% % %     A = cells{k};
% % %     Xi{k} = A(lt);   % <-- if you prefer all edges, change to A(:) and y=ymat(:)
% % % end
% % % y = ymat(lt);
% % % 
% % % M    = [Xi{:}];
% % % mask = ~(all(M==0,2) | isnan(y));
% % % M    = M(mask,:);  y = y(mask);
% % % 
% % % keep = std(M,0,1,'omitnan') > 1e-12;
% % % M    = M(:, keep);
% % % 
% % % bw        = sum( keep(1:numel(BASE_and_RED)) );  % surviving Base+RED width
% % % Xb        = M(:, 1:bw);
% % % Xf        = M;
% % % y_vec     = y;
% % % spec_cols = (bw+1):size(M,2);
% % % end
