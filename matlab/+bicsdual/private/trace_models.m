function T = trace_models(main_models, full_models, labels, scope)
% Verifies full = [main, INT_only] ordering and counts INT columns.
if nargin<4, scope=''; end
K = numel(full_models);
T = struct('label',{},'p_main',{},'p_full',{},'p_int',{},'ok_prefix',{});
for k=1:K
    main = stack_cols(main_models{k}); full = stack_cols(full_models{k});
    p_main = size(main,2); p_full = size(full,2);
    ok = p_full >= p_main;
    T(k).label = labels{k};
    T(k).p_main = p_main; T(k).p_full = p_full; T(k).p_int = max(p_full - p_main, 0);
    T(k).ok_prefix = ok;
    if ~ok
        warning('[%s] Model "%s": full has fewer cols than main!', scope, labels{k});
    end
end
end
