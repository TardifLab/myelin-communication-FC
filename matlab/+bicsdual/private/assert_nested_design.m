function dbg = assert_nested_design(Xb, Xf, spec_cols, tag, tol)
% Asserts Xf = [Xb, SPEC] and rank(Xf)>=rank(Xb). Returns a struct summary.
if nargin<5, tol = 1e-10; end
dbg = struct('ok_same_prefix',false,'rb',NaN,'rf',NaN,'dr',NaN,'p_b',NaN,'p_f',NaN,'n',size(Xb,1));
% shape checks
[pb, pf] = deal(size(Xb,2), size(Xf,2));
dbg.p_b = pb; dbg.p_f = pf;
% prefix match
ok = all(abs(Xf(:,1:pb) - Xb) < tol,'all');  % cheap but effective
dbg.ok_same_prefix = ok;
% rank check (with intercept)
rb = rank([ones(dbg.n,1) Xb]); rf = rank([ones(dbg.n,1) Xf]);
dbg.rb = rb; dbg.rf = rf; dbg.dr = rf - rb;
% spec_cols sanity
dbg.spec_cols = spec_cols(:).';
dbg.tag = tag;
if ~ok
    warning('[%s] Xf prefix != Xb (spec_cols likely misaligned).', tag);
end
if rf < rb
    warning('[%s] rank dropped (rf<rb). Design likely singular.', tag);
end
end
