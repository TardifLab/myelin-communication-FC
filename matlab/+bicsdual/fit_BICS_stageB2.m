function [SI_block, elasticity, pert] = fit_BICS_stageB2(FCs, base, tests, pinfo, Delta_block_true, opts)
% FIT_BICS_STAGEB2
% Sensitivity testing via (A) within-block edge randomization and (B) small
% multiplicative scaling of test predictors. Returns:
%   SI_block   : blocks x models x FCs sensitivity index = (ΔR²_true - mean_perm)/std_perm
%   elasticity : blocks x models x FCs ≈ d(ΔR²)/d(log scale) estimated from ±eps scaling
%   pert       : struct with summarized permutation means/stds and bookkeeping
%
% NOTES:
%     - If tests.interaction ~= 'none' and interact_at includes the level, 
%       we compare (BASE + mains) vs (BASE + mains + INT-only) with mains 
%       passed in opts.reduced_models to enforce strictly nested tests.
%     - SI is computed against the Stage A ΔR² for the exact same (block, model, band). 
%       Ensure you pass the Stage A cube from the matching interaction setting.
%
% Inputs
%   ... same as Stage A: FC cell, base & tests structs, pinfo
%   Delta_block_true : B x M x Nfc cube from Stage A (for reference)
%   opts:
%     .nperm (default 200)
%     .eps_scale (default 0.05)                               % elasticity
%     .rng_seed (optional)
%     .use_lower_only (default true) — must match Stage A
%     .perm_mode (default yresid)
%     .elasticity_mode ('zrow' (default), 'zadd', 'xscale')
%
% 2025 Mark C Nelson (MNI)
% -------------------------------------------------------------------------

if nargin<6, opts=struct; end
nperm           = getopt(opts,'nperm',200);
eps_s           = getopt(opts,'eps_scale',0.05);
use_lower       = getopt(opts,'use_lower_only',true);
seed0           = getopt(opts,'rng_seed',13);
perm_mode       = getopt(opts,'perm_mode','yresid');                        % 'Xjoint' | 'Yresid' | 'Xspec'
elasticity_mode = getopt(opts,'elasticity_mode','zrow');                    % 'zrow' (row-wise), 'zadd', 'xscale'
Kdir            = getopt(opts,'elasticity_dirs',16);                        % number of random directions (zrow)


% Unpack base
if isfield(base,'Xcells'), BASE = base.Xcells;
else, BASE = {base.Xb, base.Xc, base.Xd}; end

[models, ~] = models_from_tests(tests);
has_reduced = isfield(opts,'reduced_models') && ~isempty(opts.reduced_models);
if has_reduced
    REDUCED = opts.reduced_models; 
    assert(numel(REDUCED)==numel(models), 'reduced_models must match test models');
else
    REDUCED = cell(size(models));
end

Nfc   = numel(FCs);
Nntwk = numel(pinfo.clabels_short);

% Block map
[blk_ij, ~] = make_block_map(Nntwk, use_lower, pinfo);
assert(size(blk_ij,1) == size(Delta_block_true,1), 'Block count mismatch.');

B = size(blk_ij,1);
M = numel(models);

perm_mean = zeros(B,M,Nfc);
perm_std  = zeros(B,M,Nfc);
SI_block  = zeros(B,M,Nfc);
elasticity= zeros(B,M,Nfc);

% if isempty(gcp('nocreate')), parpool; end
% parfor f = 1:Nfc
for f = 1:Nfc 
    PM  = zeros(B,M); PS = zeros(B,M); SI = zeros(B,M); EL = zeros(B,M);
    ymat = FCs{f};
    for b = 1:B
        i = blk_ij(b,1); j = blk_ij(b,2);
        for m = 1:M
            RED = REDUCED{m};
            SPEC = models{m};

            % Collect block pack with RED included in base side
            [Xb_blk, Xf_blk, y_blk, spec_cols] = pull_block_pack2_with_spec([BASE RED], SPEC, ymat, pinfo.cis{i}, pinfo.cis{j});
            [Xb_blk, Xf_blk, y_blk] = sanitize_naninf_rows(Xb_blk, Xf_blk, y_blk);

         % ------------------------------------------------------------------
         % % % % debug (toggle off when done)
         % % %    [~,rb] = qr([ones(size(y_blk)), Xb_blk],0); rb = rank(Xb_blk);
         % % %    [~,rf] = qr([ones(size(y_blk)), Xf_blk],0); rf = rank(Xf_blk);
         % % %    fprintf('rf-rb=%d | spec_cols=%d | n=%d | p_b=%d | p_f=%d\n', rf-rb, numel(spec_cols), size(Xb_blk,1), size(Xb_blk,2), size(Xf_blk,2));
         % % % % You should always see rf - rb >= 1 when SPEC has at least one independent column.
         % ------------------------------------------------------------------

         % ------------------------------------------------------------------
         % Optional debug call
            if isfield(opts,'debug') && isfield(opts.debug,'probe_B') && opts.debug.probe_B
                probe_compare_snapshot(y_blk, Xb_blk, Xf_blk, spec_cols, [], [], ...
                   sprintf('B2-BLOCK b=%d m=%d f=%d', b, m, f), getopt(opts.debug,'save',''));
            end
          % ------------------------------------------------------------------

            if isempty(y_blk) || numel(y_blk) < (size(Xf_blk,2) + 5)
                PM(b,m)=0; PS(b,m)=0; SI(b,m)=0; EL(b,m)=0; continue;
            end

            % Permutations
            rng(seed0 + f*1e6 + b*1e3 + m, 'twister');
            null_vals = zeros(nperm,1);
            switch lower(perm_mode)
                case 'xjoint'
                    for p=1:nperm
                        idx = randperm(numel(y_blk));
                        Xb_perm = Xb_blk(idx,:); Xf_perm = Xf_blk(idx,:);
                        dR2 = safe_compare3(Xb_perm, Xf_perm, y_blk); if ~isfinite(dR2), dR2=0; end
                        null_vals(p) = dR2;
                    end
                case 'xspec'
                    for p=1:nperm
                        idx = randperm(numel(y_blk));
                        Xf_perm = Xf_blk; Xf_perm(:,spec_cols)=Xf_perm(idx,spec_cols);
                        dR2 = safe_compare3(Xb_blk, Xf_perm, y_blk); if ~isfinite(dR2), dR2=0; end
                        null_vals(p) = dR2;
                    end
                case 'yresid'
                    % Mred   = fitlm(Xb_blk, y_blk);
                    % yhat_r = Mred.Fitted; eres = Mred.Residuals.Raw;
                    [yhat_r, eres] = reduced_fit_qr(y_blk, Xb_blk);         % QR version
                    for p=1:nperm
                        idx = randperm(numel(eres));
                        y_perm = yhat_r + eres(idx);
                        dR2 = safe_compare3(Xb_blk, Xf_blk, y_perm); if ~isfinite(dR2), dR2=0; end
                        null_vals(p) = dR2;
                    end
                otherwise
                    error('Unknown perm_mode %s', perm_mode);
            end

            mu = mean(null_vals); sd = std(null_vals) + eps;
            PM(b,m)=mu; PS(b,m)=sd;
            d_true = Delta_block_true(b,m,f);
            SI(b,m) = (d_true - mu) / sd;

        % ----- (B) Elasticity on SPEC only --------------------------------------
        % Modes:
        %   'zrow'  (best for lin reg): row-wise heterogenous dilation of SPEC columns.
        %                          For each of K random directions w (per column),
        %                          X_up = X .* (1 + eps_s * w), X_dn = X .* (1 - eps_s * w)
        %                          denom = 2*eps_s. Returns the mean over K directions.
        %   'zadd'  (center-preserving column scaling): X_up = (1+eps_s)*X, X_dn = (1-eps_s)*X
        %            -> ΔR^2 invariant to per-column scaling => ~0 (kept for compatibility).
        %   'xscale' (multiplicative): same scaling as 'zadd' but per-log denominator.
            if isempty(spec_cols)
                EL(b,m) = 0;
            else
                switch lower(elasticity_mode)
            
                    case 'zrow'
                        % Deterministic RNG per (f,b,m) so it's reproducible
                        rng(seed0 + f*1e6 + b*1e3 + m, 'twister');
            
                        % Precompute column-wise row counts
                        nobs = size(Xf_blk,1);
            
                        dir_vals = zeros(Kdir,1);
                        for kk = 1:Kdir
                            Xup = Xf_blk; Xdn = Xf_blk;
            
                            % Build a row-wise, zero-mean, unit-sd jitter per SPEC column
                            % Each column gets its own independent direction vector w
                            for c = spec_cols
                                w = randn(nobs,1);
                                w = (w - mean(w)) / max(std(w,0,1), 1e-12);   % standardize
                                Xup(:,c) = Xup(:,c) .* (1 + eps_s * w);
                                Xdn(:,c) = Xdn(:,c) .* (1 - eps_s * w);
                            end
            
                            d_up = safe_compare3(Xb_blk, Xup, y_blk); if ~isfinite(d_up), d_up = 0; end
                            d_dn = safe_compare3(Xb_blk, Xdn, y_blk); if ~isfinite(d_dn), d_dn = 0; end
            
                            denom = 2 * eps_s; if denom <= 0 || ~isfinite(denom), denom = eps; end
                            dir_vals(kk) = (d_up - d_dn) / denom;
                        end
            
                        % Average directional derivatives (you could also take median)
                        EL(b,m) = mean(dir_vals, 'omitnan');
            
                    case 'zadd'
                        % Column scaling -> ΔR^2 invariant; but could be useful in other contexts?
                        Xup = Xf_blk; Xdn = Xf_blk;
                        Xup(:,spec_cols) = (1+eps_s) * Xup(:,spec_cols);
                        Xdn(:,spec_cols) = (1-eps_s) * Xdn(:,spec_cols);
            
                        d_up = safe_compare3(Xb_blk, Xup, y_blk); if ~isfinite(d_up), d_up = 0; end
                        d_dn = safe_compare3(Xb_blk, Xdn, y_blk); if ~isfinite(d_dn), d_dn = 0; end
            
                        denom = 2 * eps_s; if denom <= 0 || ~isfinite(denom), denom = eps; end
                        EL(b,m) = (d_up - d_dn) / denom;
            
                    case 'xscale'
                        % Multiplicative + log denominator (if X not z-scored)
                        Xup = Xf_blk; Xdn = Xf_blk;
                        Xup(:,spec_cols) = (1+eps_s) * Xup(:,spec_cols);
                        Xdn(:,spec_cols) = (1-eps_s) * Xdn(:,spec_cols);
            
                        d_up = safe_compare3(Xb_blk, Xup, y_blk); if ~isfinite(d_up), d_up = 0; end
                        d_dn = safe_compare3(Xb_blk, Xdn, y_blk); if ~isfinite(d_dn), d_dn = 0; end

                        denom = log(1+eps_s) - log(1-eps_s);
                        if denom <= 0 || ~isfinite(denom), denom = eps; end
                        EL(b,m) = (d_up - d_dn) / denom;
            
                    otherwise
                        error('Unknown opts.elasticity_mode: %s (use ''zrow'', ''zadd'' or ''xscale'')', elasticity_mode);
                end
            end
        end
    end
    perm_mean(:,:,f)=PM; perm_std(:,:,f)=PS; SI_block(:,:,f)=SI; elasticity(:,:,f)=EL;
end

pert = struct('perm_mean',perm_mean,'perm_std',perm_std,'nperm',nperm,'eps_scale',eps_s);
% -------------------------------------------------------------------------
end
