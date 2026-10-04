function [Delta_all, P_all, stats_all] = fit_BICS_dual_stageA(ALLY, ALLX, CC, tests, pinfo, opts)
% FIT_BICS_DUAL_STAGEA
% Orchestrates Level-1 (Route/Diff), Level-2 (Both)
% Interactions can be included with caliber or ED
% using generalized Stage A that supports arbitrary base width.
%
% Inputs
%   ALLY   : 1xNfc cell of NxN FC matrices
%   ALLX   : 1x6 cell (order: {caliber, MTsat, gratio, delay, binary, ED}), each 1x6 models, except ED N x N
%   CC     : [iRoute iDiff] (indices 1..6), e.g., [1 6] = SPE × DE
%   tests  : struct with fields:
%            - myelin_predictors    : cellstr subset of {'MTsat','gratio','delay'}
%            - modes                : cellstr subset of {'single','pairs','all'}
%            - levels               : e.g., [1 2]
%            - label_comm_models    : e.g., 'SPE-CMY'
%   pinfo  : struct with parcel info e.g., .cis, .clabels_short, etc.
%   opts   : struct (optional):
%            - alpha (0.05)
%            - use_lower_only (true)
%
% Outputs
%   Delta_all : struct with fields per level:
%       .L1_route, .L1_diff, .L2_both (each B x M x Nfc)
%   P_all     : same structure as Delta_all for p-values
%   stats_all : struct with metadata for each level and merged labels
%
% 2025 Mark C Nelson (MNI)
% -------------------------------------------------------------------------


if nargin<6, optsA = struct; else, optsA = opts; end
% label_model = getopt(opts,'label_model','');

% Validate CC (one routing-class, one diffusion-class)
routing_mask  = [1 1 1 1 0 0];  % SPE NE SIE PT CMY DE
diffusion_mask= [0 0 0 0 1 1];
iR = CC(1); iD = CC(2);
assert( (routing_mask(iR)==1 && diffusion_mask(iD)==1) || ...
        (routing_mask(iD)==1 && diffusion_mask(iR)==1), ...
        'CC must contain one routing-class and one diffusion-class model.');
% normalize ordering: R first, D second
if diffusion_mask(iR)==1, tmp=iR; iR=iD; iD=tmp; end

% Resolve predictor matrices from useCM
U = unpack_useCM(ALLX);   % struct with fields: caliber{6}, MTsat{6}, gratio{6}, delay{6}, binary{6}

% Build base blocks
ED = U.ED; 
if isempty(ED)
    error('Please set U.ED inside unpack_useCM() or adapt to pass ED explicitly.');
end

% Hardcoded original version
% Base_route_X = { U.binary{iR}, U.caliber{iR}, ED };       % 3 base cols
% Base_diff_X  = { U.binary{iD}, U.caliber{iD}, ED };
% Base_both_X  = { U.binary{iR}, U.caliber{iR}, ED, U.binary{iD}, U.caliber{iD} }; % 5 base cols

% Build myelin menus per CM based on tests.myelin_predictors
mp                  = getopt(tests,'myelin_predictors',{'MTsat','gratio','delay'});
base_predictors     = canonical_base_predictors(tests);
modes               = getopt(tests,'modes',{'single','pairs','all'});
levels              = getopt(tests,'levels',[1 2]);
comm_model_labels   = getopt(tests,'label_comm_models','R-D');
interaction         = getopt(tests, 'interaction', 'none');                 % 'none'|'caliber'|'ed'|'both'
interact_at         = getopt(tests, 'interact_at', 'both');                 % 'L1'|'L2'|'both'
effect_mode         = getopt(tests, 'effect_mode', 'incremental');          % 'total'|'incremental'

% Gather R and D base stacks
[Base_route_X, Base_route_labels] = base_pack(U, iR, base_predictors, 'R');
[Base_diff_X,  Base_diff_labels]  = base_pack(U, iD, base_predictors, 'D');
Base_both_X         = [Base_route_X  {ED}  Base_diff_X];                    % Add ED to all base packs  
Base_route_X        = [Base_route_X, {ED}];
Base_diff_X         = [Base_diff_X, {ED}];
Base_both_labels    = [Base_route_labels 'ED' Base_diff_labels];            % adjust labels for ED
Base_route_labels   = [Base_route_labels 'ED'];                             
Base_diff_labels    = [Base_diff_labels  'ED'];

% Gather R and D myelin stacks in a canonical order matching labels
% [R_Xlist, R_labels] = myelin_pack(U, iR, mp);
% [D_Xlist, D_labels] = myelin_pack(U, iD, mp);
[R_Xlist, R_labels] = predictor_pack(U, iR, mp);
[D_Xlist, D_labels] = predictor_pack(U, iD, mp);

% Build model menus per modes
[R_models, R_mlabels] = build_models_strict(R_Xlist, R_labels, modes);
[D_models, D_mlabels] = build_models_strict(D_Xlist, D_labels, modes);

% ---------------------------------------------------
% Optional debugging
if isfield(optsA,'debug') && getopt(optsA.debug,'trace_menus',false)
    TR_R  = trace_models(R_models,  R_models,  R_mlabels,  'L1-R mains'); %#ok<NASGU>
    TR_D  = trace_models(D_models,  D_models,  D_mlabels,  'L1-D mains'); %#ok<NASGU>
    [RD_models, RD_labels] = cartesian_models(R_models,R_mlabels,D_models,D_mlabels);
    TR_RD = trace_models(RD_models, RD_models, RD_labels, 'L2 mains'); %#ok<NASGU>
    % If interactions requested:
    if ~strcmpi(interaction,'none')
        switch lower(interact_at)
            case 'l1'
                R_INT  = apply_interactions_flagged(R_models,'R',interaction,Base_route_X);
                D_INT  = apply_interactions_flagged(D_models,'D',interaction,Base_diff_X);
                TR_RI  = trace_models(R_models, R_INT, R_mlabels, 'L1-R full'); %#ok<NASGU>
                TR_DI  = trace_models(D_models, D_INT, D_mlabels, 'L1-D full'); %#ok<NASGU>
            case 'l2'
                RD_INT = apply_interactions_flagged(RD_models,'RD',interaction,Base_both_X);
                TR_RDI = trace_models(RD_models, RD_INT, RD_labels, 'L2 full'); %#ok<NASGU>
            case 'both'    
                R_INT  = apply_interactions_flagged(R_models,'R',interaction,Base_route_X);
                D_INT  = apply_interactions_flagged(D_models,'D',interaction,Base_diff_X);
                TR_RI  = trace_models(R_models, R_INT, R_mlabels, 'L1-R full'); %#ok<NASGU>
                TR_DI  = trace_models(D_models, D_INT, D_mlabels, 'L1-D full'); %#ok<NASGU>
                RD_INT = apply_interactions_flagged(RD_models,'RD',interaction,Base_both_X);
                TR_RDI = trace_models(RD_models, RD_INT, RD_labels, 'L2 full'); %#ok<NASGU>
        end
    end
end
% ----------------------------------------------------------


% ----- Level-1
Delta_all = struct(); P_all = struct(); stats_all = struct();
if any(levels==1)

  % --- ROUTING ---
    base_route = struct('Xcells', {Base_route_X}, 'labels_base', {Base_route_labels});
    tests_R    = tests_from_models(R_models, R_mlabels);
    o = overlay_opts(optsA, struct('label_model', 'L1-Route'));
    if strcmpi(effect_mode,'incremental')                  % incremental vs total
      % incremental gain path
        if strcmpi(interaction,'none') || ismember(lower(interact_at), {'l2'})  % no ints at L1
            [D1r, P1r, S1r] = bicsdual.fit_BICS_stageA2(ALLY, base_route, tests_R, pinfo, o);
        else
            R_INT = apply_interactions_flagged(R_models, 'R', interaction, Base_route_X);
            INT_only = cellfun(@(main, full) full(numel(main)+1:end),R_models, R_INT, 'UniformOutput', false);
            tests_L1r = tests_from_models(INT_only, R_mlabels);
            o.reduced_models = R_models; 
            o.reduced_labels = R_mlabels; 
            o.label_model    = sprintf('L1-Route+Int{%s}', lower(interaction));
            [D1r, P1r, S1r] = bicsdual.fit_BICS_stageA2(ALLY, base_route, tests_L1r, pinfo, o);
        end  
        label_out = o.label_model;
    else
      % Total contribution path
        switch lower(interaction)
            case 'none'
                % SPEC = myelin mains only
                SPEC_cells = R_models;     % each is the main-effects column set (beyond base)

            case {'caliber','ed','both'}
                % Build interactions variants to extract the *added* interaction columns.
                do_cal = ismember(lower(interaction), {'caliber','both'}); %any(strcmpi(interaction, {'caliber','both'}));
                do_ed  = ismember(lower(interaction), {'ed','both'});

                % Start from mains:
                SPEC_cells = R_models;

                % CALIBER interactions (if requested at L1/both)
                if do_cal && (strcmpi(interact_at,'both') || strcmpi(interact_at,'l1'))
                    R_INT_cal = apply_interactions_flagged(R_models, 'R', 'caliber', Base_route_X);
                    INT_cal_only = cellfun(@(main,full) full(numel(main)+1:end), R_models, R_INT_cal, 'UniformOutput', false);
                    SPEC_cells = cellfun(@(spec,add) [spec, add], SPEC_cells, INT_cal_only, 'UniformOutput', false);
                end

                % ED interactions (if requested at L1/both)
                if do_ed && (strcmpi(interact_at,'both') || strcmpi(interact_at,'l1'))
                    R_INT_ed = apply_interactions_flagged(R_models, 'R', 'ed', Base_route_X);
                    INT_ed_only = cellfun(@(main,full) full(numel(main)+1:end), R_models, R_INT_ed, 'UniformOutput', false);
                    SPEC_cells = cellfun(@(spec,add) [spec, add], SPEC_cells, INT_ed_only, 'UniformOutput', false);
                end

            otherwise
                error('Unknown interaction option: %s', interaction);
        end

        % Turn SPEC_cells into a tests struct (each entry = what we add on top of base)
        tests_total = tests_from_models(SPEC_cells, R_mlabels);

        % Label and run (no reduced_models => RED=base; FULL=base+SPEC)
        o2 = o;
        o2 = overlay_opts(o2, struct('label_model', sprintf('L1-Route TOTAL{%s}', lower(interaction))));
        % IMPORTANT: DO NOT pass o2.reduced_models here
        if isfield(o2,'reduced_models'); o2 = rmfield(o2, {'reduced_models','reduced_labels'}); end
        [D1r, P1r, S1r] = bicsdual.fit_BICS_stageA2(ALLY, base_route, tests_total, pinfo, o2);
        label_out = o2.label_model;
    end

    Delta_all.L1_route = D1r; P_all.L1_route = P1r; stats_all.L1_route = annotate_stats(S1r, label_out);


  % --- DIFFUSION ---
    base_diff  = struct('Xcells', {Base_diff_X}, 'labels_base', {Base_diff_labels});
    tests_D    = tests_from_models(D_models, D_mlabels);
    o = overlay_opts(optsA, struct('label_model', 'L1-Diff'));
    if strcmpi(effect_mode,'incremental')                  % incremental vs total
      % incremental gain path
        if strcmpi(interaction,'none') || ismember(lower(interact_at), {'l2'})
            [D1d, P1d, S1d] = bicsdual.fit_BICS_stageA2(ALLY, base_diff, tests_D, pinfo, o);
        else
            D_INT = apply_interactions_flagged(D_models, 'D', interaction, Base_diff_X);
            INT_only = cellfun(@(main, full) full(numel(main)+1:end),D_models,D_INT,'UniformOutput',false);
            tests_L1d = tests_from_models(INT_only, D_mlabels);
            o.reduced_models = D_models; 
            o.reduced_labels = D_mlabels; 
            o.label_model    = sprintf('L1-Diff+Int{%s}', lower(interaction));
            [D1d, P1d, S1d] = bicsdual.fit_BICS_stageA2(ALLY, base_diff, tests_L1d, pinfo, o);
        end
        label_out = o.label_model;
    else
      % Total contribution path
        switch lower(interaction)
            case 'none'
                % SPEC = myelin mains only
                SPEC_cells = D_models;     % each is the main-effects column set (beyond base)

            case {'caliber','ed','both'}
                % Build interactions variants to extract the *added* interaction columns.
                do_cal = ismember(lower(interaction), {'caliber','both'}); %any(strcmpi(interaction, {'caliber','both'}));
                do_ed  = ismember(lower(interaction), {'ed','both'});

                % Start from mains:
                SPEC_cells = D_models;

                % CALIBER interactions (if requested at L1/both)
                if do_cal && (strcmpi(interact_at,'both') || strcmpi(interact_at,'l1'))
                    D_INT_cal = apply_interactions_flagged(D_models, 'D', 'caliber', Base_diff_X);
                    INT_cal_only = cellfun(@(main,full) full(numel(main)+1:end), D_models, D_INT_cal, 'UniformOutput', false);
                    SPEC_cells = cellfun(@(spec,add) [spec, add], SPEC_cells, INT_cal_only, 'UniformOutput', false);
                end

                % ED interactions (if requested at L1/both)
                if do_ed && (strcmpi(interact_at,'both') || strcmpi(interact_at,'l1'))
                    D_INT_ed = apply_interactions_flagged(D_models, 'D', 'ed', Base_diff_X);
                    INT_ed_only = cellfun(@(main,full) full(numel(main)+1:end), D_models, D_INT_ed, 'UniformOutput', false);
                    SPEC_cells = cellfun(@(spec,add) [spec, add], SPEC_cells, INT_ed_only, 'UniformOutput', false);
                end

            otherwise
                error('Unknown interaction option: %s', interaction);
        end

        % Turn SPEC_cells into a tests struct (each entry = what we add on top of base)
        tests_total = tests_from_models(SPEC_cells, D_mlabels);

        % Label and run (no reduced_models => RED=base; FULL=base+SPEC)
        o2 = o;
        o2 = overlay_opts(o2, struct('label_model', sprintf('L1-Diff TOTAL{%s}', lower(interaction))));
        % IMPORTANT: DO NOT pass o2.reduced_models here
        if isfield(o2,'reduced_models'); o2 = rmfield(o2, {'reduced_models','reduced_labels'}); end
        [D1d, P1d, S1d] = bicsdual.fit_BICS_stageA2(ALLY, base_diff, tests_total, pinfo, o2);
        label_out = o2.label_model;
    end
    Delta_all.L1_diff = D1d; P_all.L1_diff = P1d; stats_all.L1_diff = annotate_stats(S1d, label_out);
end


% ----- Level-2 BOTH (Cartesian of R and D model menus)
if any(levels>=2)
    [RD_models, RD_labels] = cartesian_models(R_models, R_mlabels, D_models, D_mlabels);
    base_both = struct('Xcells', {Base_both_X}, 'labels_base', {Base_both_labels});
    tests_RD  = tests_from_models(RD_models, RD_labels);
    o = overlay_opts(optsA, struct('label_model', 'L2-Both'));
    if strcmpi(effect_mode,'incremental')                  % incremental vs total
      % incremental gain path
        if strcmpi(interaction,'none') || ismember(lower(interact_at), {'l1'})  % no new ints at L2
            [D2, P2, S2] = bicsdual.fit_BICS_stageA2(ALLY, base_both, tests_RD, pinfo, o);
        else
            % Apply interactions to BOTH halves; reduced = RD main-effects model
            RD_INT = apply_interactions_flagged(RD_models, 'RD', interaction, Base_both_X);
            INT_only = cellfun(@(main,full) full(numel(main)+1:end),RD_models, RD_INT,'UniformOutput',false);
            tests_L2 = tests_from_models(INT_only, RD_labels);
            o.reduced_models = RD_models; 
            o.reduced_labels = RD_labels; 
            o.label_model = sprintf('L2-Both+Int{%s}', lower(interaction));
            [D2, P2, S2] = bicsdual.fit_BICS_stageA2(ALLY, base_both, tests_L2, pinfo, o);
        end
    else
      % Total contribution path
        switch lower(interaction)
            case 'none'
                % SPEC = myelin mains only
                SPEC_cells = RD_models;     % each is the main-effects column set (beyond base)

            case {'caliber','ed','both'}
                % Build interactions variants to extract the *added* interaction columns.
                do_cal = any(strcmpi(interaction, {'caliber','both'}));
                do_ed  = any(strcmpi(interaction, {'ed','both'}));

                % Start from mains:
                SPEC_cells = RD_models;

                % CALIBER interactions (if requested at L2/both)
                if do_cal && (strcmpi(interact_at,'both') || strcmpi(interact_at,'l2'))
                    RD_INT_cal = apply_interactions_flagged(RD_models, 'RD', 'caliber', Base_both_X);
                    INT_cal_only = cellfun(@(main,full) full(numel(main)+1:end), RD_models, RD_INT_cal, 'UniformOutput', false);
                    SPEC_cells = cellfun(@(spec,add) [spec, add], SPEC_cells, INT_cal_only, 'UniformOutput', false);
                end

                % ED interactions (if requested at L2/both)
                if do_ed && (strcmpi(interact_at,'both') || strcmpi(interact_at,'l2'))
                    RD_INT_ed = apply_interactions_flagged(RD_models, 'RD', 'ed', Base_both_X);
                    INT_ed_only = cellfun(@(main,full) full(numel(main)+1:end), RD_models, RD_INT_ed, 'UniformOutput', false);
                    SPEC_cells = cellfun(@(spec,add) [spec, add], SPEC_cells, INT_ed_only, 'UniformOutput', false);
                end

            otherwise
                error('Unknown interaction option: %s', interaction);
        end

        % Turn SPEC_cells into a tests struct (each entry = what we add on top of base)
        tests_total = tests_from_models(SPEC_cells, RD_labels);

        % Label and run (no reduced_models => RED=base; FULL=base+SPEC)
        o2 = o;
        o2 = overlay_opts(o2, struct('label_model', sprintf('L2-Both TOTAL{%s}', lower(interaction))));
        % IMPORTANT: ensure we do NOT pass o2.reduced_models here
        if isfield(o2,'reduced_models'); o2 = rmfield(o2, {'reduced_models','reduced_labels'}); end

        [D2, P2, S2] = bicsdual.fit_BICS_stageA2(ALLY, base_both, tests_total, pinfo, o2);
    end

    Delta_all.L2_both = D2; P_all.L2_both = P2; stats_all.L2_both = annotate_stats(S2, o.label_model);
end

% % % % Level-3 Interactions (optional)
% % % if any(levels==3) && doInt
% % %     assert(isfield(Delta_all,'L2_both'), 'L3 requires Level-2 results.');
% % %     % Build L3 by adding interactions to each RD model in the SAME order
% % %     [RD_models, RD_labels] = cartesian_models(R_models, R_mlabels, D_models, D_mlabels);
% % %     [INT_models, INT_labels] = add_interactions_on_RDboth(RD_models, RD_labels, Base_both_X, {'caliber_R','caliber_D'}); %#ok<ASGLU>
% % %     % For L3, reduced model = Base_both + (R main + D main) = the corresponding L2 model
% % %     tests_L3 = tests_from_models(INT_models, INT_labels);
% % %     opts.reduced_models = RD_models;       opts.reduced_labels = RD_labels;
% % %     o = overlay_opts(optsA, struct('base_cols', numel(Base_both_X), 'label_model', 'L3-Both+Int', ...
% % %                                    'reduced_models', RD_models, 'reduced_labels', RD_labels ));
% % %     [D3, P3, S3] = bicsdual.fit_BICS_stageA2(ALLY, base_both, tests_L3, pinfo, o);
% % %     Delta_all.L3_int = D3; P_all.L3_int = P3; stats_all.L3_int = annotate_stats(S3, 'L3-Both+Int');
% % % end

% Attach top-level label
stats_all.meta = struct('label_comm_models', comm_model_labels, ...
                        'CC',[iR iD], ...
			'base_predictors',{base_predictors}, ...
                        'effect_mode',effect_mode, ...
                        'interaction_label', ['int-' interaction]);
% -------------------------------------------------------------------------
end
