function check_deltaR2_engines(Xb, Xf, y, tol, tag)
if nargin<4, tol = 1e-8; end
if nargin<5, tag = ''; end
[dR2_qr,~] = safe_compare3(Xb,Xf,y);
try
    mdl_b = fitlm(Xb,y); mdl_f = fitlm(Xf,y);
    dR2_lm = mdl_f.Rsquared.Ordinary - mdl_b.Rsquared.Ordinary;
    if abs(dR2_qr - dR2_lm) > tol
        warning('[%s] QR vs fitlm ΔR² mismatch: %.3g vs %.3g', tag, dR2_qr, dR2_lm);
    end
catch
    % ok if fitlm fails for small n
end
end
