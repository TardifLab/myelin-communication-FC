function [Delta_block, P_block, stats] = fit_BICS_stageA2(FCs, base, tests, pinfo, opts)
% FIT_BICS_STAGEA2
% Computes ΔR² per RSN block, model spec, and FC target.
%
% NOTES: 
%     - If opts.reduced_models is supplied, they are treated as "mains" (RED),
%       and SPEC is the INT-only part. Collectors *_with_spec concatenate [BASE RED]
%       on the left and SPEC on the right so compare is (BASE+RED) vs (BASE+RED+SPEC).
%     - If tests.interaction ~= 'none' and interact_at includes the level, 
%       we compare (BASE + mains) vs (BASE + mains + INT-only) with mains 
%       passed in opts.reduced_models to enforce strictly nested tests.
%
%
% Inputs
%   FCs    : 1xNfc cell of NxN FC matrices 
%   base   : struct with fields:
%             - accepts base.Xcells (cell array) of arbitrary length
%   tests  : struct with fields:
%            - X1 (NxN), X2 (NxN), X3 (NxN) : MTsat, g-ratio, delay
%            - labels_tests : 1x3 cellstr, e.g., {'MTsat','g-ratio','delay'}
%   pinfo  : struct with fields 'clabels_short' & 'cis'
%   opts   : struct (optional):
%            - alpha (0.05)
%            - use_lower_only (true)
%            - label_mode : model level and data type (e.g., L1-Route)
%            - optionally accepts opts.reduced_models to compare (Base + ReducedSpec) vs (Base + ReducedSpec + TestSpec)
%
% Outputs
%   Delta_block : B x M x Nfc
%   P_block     : B x M x Nfc
%   stats       : struct with global/ntwk/node cubes and labels
%
% 2025 Mark C Nelson (MNI)
% -------------------------------------------------------------------------


if nargin<5, opts=struct; end
alpha           = getopt(opts,'alpha',0.05);
use_lower       = getopt(opts,'use_lower_only',true);
label_model     = getopt(opts,'label_model','');
% base_cols       = getopt(opts,'base_cols',3);
hasRED          = isfield(opts,'reduced_models') && ~isempty(opts.reduced_models);

% Unpack base
if isfield(base,'Xcells'), BASE = base.Xcells;
else, BASE = {base.Xb, base.Xc, base.Xd}; end
LBL_base = ensure_labels(base,'labels_base',arrayfun(@(i) sprintf('base%02d',i),1:numel(BASE),'uni',0));

% Build tests list
[models, model_labels] = models_from_tests(tests);

% Reduced models (for INTs)
has_reduced = isfield(opts,'reduced_models') && ~isempty(opts.reduced_models);
if has_reduced
    REDUCED = opts.reduced_models;  % cell array (length = M), each reduced SPEC (cell of matrices)
    assert(numel(REDUCED)==numel(models), 'reduced_models must match number of test models');
else
    REDUCED = cell(size(models)); % all empty
end

% Dimensions
Nfc   = numel(FCs);
Nnode = size(FCs{1},1);
Nntwk = numel(pinfo.clabels_short);
M     = numel(models);

% Block map
[blk_ij, blk_labels] = make_block_map(Nntwk, use_lower, pinfo);
B = size(blk_ij,1);

% Prealloc
Delta_block = zeros(B, M, Nfc);
P_block     = ones (B, M, Nfc);
global_deltaR2 = zeros(M, Nfc); global_p = ones(M, Nfc);
ntwk_deltaR2   = zeros(M, Nntwk, Nntwk, Nfc); ntwk_p = ones(M, Nntwk, Nntwk, Nfc);
node_deltaR2   = zeros(M, Nnode, Nfc);        node_p = ones(M, Nnode, Nfc);

% ===== GLOBAL + NODE (serial) =====
for m = 1:M
    SPEC = models{m};
    RED  = REDUCED{m};
    for f = 1:Nfc
        ymat = FCs{f};
        [Xb,Xf,yv] = collect_edges2_with_spec([BASE RED], SPEC, ymat);
         % ------------------------------------------------------------------
         % Optional debug call: V1
            if isfield(opts,'debug') && isfield(opts.debug,'probe') && opts.debug.probe
                probe_compare_snapshot(yv, Xb, Xf, size(Xb,2)+(1:(size(Xf,2)-size(Xb,2))), ...
                    annotate_names(BASE, RED), annotate_names(base, { [RED(:).', SPEC(:).'] }), ...
                    sprintf('A2-GLOBAL m=%d f=%d', m, f), getopt(opts.debug,'save',''));
            end
         % ------------------------------------------------------------------

        [Xb, Xf, yv] = sanitize_naninf_rows(Xb, Xf, yv);                    % dim=1 default
        [dR2, pval] = safe_compare3(Xb, Xf, yv);
        global_deltaR2(m,f) = dR2; global_p(m,f) = pval;

        for n = 1:Nnode
            [Xb_n, Xf_n, y_n] = pull_node_pack2_with_spec([BASE RED], SPEC, ymat, n);
            [Xb_n, Xf_n, y_n] = sanitize_naninf_rows(Xb_n, Xf_n, y_n);
            [dR2, pval] = safe_compare3(Xb_n, Xf_n, y_n);
            node_deltaR2(m,n,f) = dR2; node_p(m,n,f) = pval;
        end
    end
end

% ===== BLOCKS (parallel over FCs) =====
% if isempty(gcp('nocreate')), parpool; end
% parfor f = 1:Nfc
for f = 1:Nfc 
    Dbf = zeros(B, M);  Pbf = ones(B, M);
    ymat = FCs{f};
    for b = 1:B
        i = blk_ij(b,1); j = blk_ij(b,2);
        io = pinfo.cis{i}; jo = pinfo.cis{j};
        for m = 1:M
            SPEC = models{m}; RED = REDUCED{m};
            [Xb_blk, Xf_blk, y_blk, spec_cols] = pull_block_pack2_with_spec([BASE RED], SPEC, ymat, io, jo);
            [Xb_blk, Xf_blk, y_blk] = sanitize_naninf_rows(Xb_blk, Xf_blk, y_blk);

         % ------------------------------------------------------------------
         % Optional debug call: V1
            if isfield(opts,'debug') && isfield(opts.debug,'probe_blocks') && opts.debug.probe_blocks
                probe_compare_snapshot(y_blk, Xb_blk, Xf_blk, size(Xb_blk,2)+(1:(size(Xf_blk,2)-size(Xb_blk,2))), ...
                    annotate_names(BASE, RED), annotate_names([BASE RED SPEC]), ...
                    sprintf('A2-BLOCK b=%d m=%d f=%d', b, m, f), getopt(opts.debug,'save',''));
            end
         % ------------------------------------------------------------------
          % Model
            [dR2, pval] = safe_compare3(Xb_blk, Xf_blk, y_blk);
            Dbf(b,m) = dR2;  Pbf(b,m) = pval;
        end
    end
    Delta_block(:,:,f) = Dbf;
    P_block   (:,:,f) = Pbf;
end

% Rebuild RSN×RSN tensors
for f = 1:Nfc
    for m = 1:M
        for b = 1:B
            i = blk_ij(b,1); j = blk_ij(b,2);
            d = Delta_block(b,m,f); p = P_block(b,m,f);
            ntwk_deltaR2(m,i,j,f) = d; ntwk_p(m,i,j,f) = p;
            ntwk_deltaR2(m,j,i,f) = d; ntwk_p(m,j,i,f) = p;
        end
    end
end

% Pack stats
stats = struct();
stats.global_deltaR2 = global_deltaR2;
stats.global_p       = global_p;
stats.ntwk_deltaR2   = ntwk_deltaR2;
stats.ntwk_p         = ntwk_p;
stats.node_deltaR2   = node_deltaR2;
stats.node_p         = node_p;
stats.model_labels   = model_labels;
stats.block_map      = blk_ij;
stats.block_labels   = blk_labels;
stats.use_lower      = use_lower;
stats.Nntwk          = Nntwk;
stats.alpha          = alpha;
if ~isempty(label_model), stats.label_model = label_model; end
% -------------------------------------------------------------------------
end
