function REP = probe_compare_snapshot(y, Xb, Xf, spec_cols, names_b, names_f, tag, save_mat)
% One-shot snapshot: asserts nestedness, prints collinearity, compares ΔR² via safe_compare3 vs fitlm (optional).
if nargin<8, save_mat=''; end
if nargin<7, tag=''; end
if nargin<6 || isempty(names_f), names_f = arrayfun(@(i) sprintf('Xf%02d',i),1:size(Xf,2),'uni',0); end
if nargin<5 || isempty(names_b), names_b = arrayfun(@(i) sprintf('Xb%02d',i),1:size(Xb,2),'uni',0); end

REP = struct();
REP.tag = tag;
REP.n = numel(y); REP.p_b = size(Xb,2); REP.p_f = size(Xf,2);
REP.spec_cols = spec_cols(:).';
REP.assert = assert_nested_design(Xb, Xf, spec_cols, tag);

% ΔR² (QR)
[REP.dR2_qr, REP.p_qr] = safe_compare3(Xb, Xf, y);

% ΔR² (fitlm cross-check — optional; guard tiny n)
try
    mdl_b = fitlm(Xb, y); mdl_f = fitlm(Xf, y);
    R2_b  = mdl_b.Rsquared.Ordinary; R2_f = mdl_f.Rsquared.Ordinary;
    REP.dR2_lm = R2_f - R2_b;
catch
    REP.dR2_lm = NaN;
end

% Design signatures
REP.sig_b = design_signature(Xb, names_b);
REP.sig_f = design_signature(Xf, names_f);

% Collinearity report (full)
REP.colrep = collinearity_report(Xf, names_f);

% Print brief line
fprintf('[%s] n=%d | p_b=%d p_f=%d | dr=%d | dR2_qr=%.4g | lm=%.4g\n', ...
        tag, REP.n, REP.p_b, REP.p_f, REP.assert.dr, REP.dR2_qr, REP.dR2_lm);

if ~isempty(save_mat)
    save(save_mat, 'REP');
end
end
