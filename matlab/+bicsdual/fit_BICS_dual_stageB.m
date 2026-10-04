function [SI_all, elasticity_all, pert_all] = fit_BICS_dual_stageB(ALLY, ALLX, CC, tests, pinfo, Delta_all, opts)
% ORCHESTRATION: STAGE B (Sensitivity / Elasticity) for dual-CM hierarchy
% Mirrors Stage A orchestration, calling StageB2 with matching reduced models.
%
% 2025 Mark C Nelson (MNI)
% -------------------------------------------------------------------------

if nargin<7, optsB = struct; else, optsB = opts; end

% Rebuild the same menus so we pass consistent reduced models to StageB2
routing_mask  = [1 1 1 1 0 0]; diffusion_mask= [0 0 0 0 1 1];
iR = CC(1); iD = CC(2);
if diffusion_mask(iR)==1, tmp=iR; iR=iD; iD=tmp; end

U = unpack_useCM(ALLX);
ED = U.ED; if isempty(ED), error('Please set U.ED in unpack_useCM().'); end
canonical_base_predictors(tests);
Base_route_X = { U.binary{iR}, U.caliber{iR}, ED };
Base_diff_X  = { U.binary{iD}, U.caliber{iD}, ED };
Base_both_X  = { U.binary{iR}, U.caliber{iR}, ED, U.binary{iD}, U.caliber{iD} };

mp                  = getopt(tests,'myelin_predictors',{'myln1','myln2','myln3'});
modes               = getopt(tests,'modes',{'single','pairs','all'});
levels              = getopt(tests,'levels',[1 2]);
interaction         = getopt(tests, 'interaction', 'none');                 % 'none'|'caliber'|'ed'|'both'
interact_at         = getopt(tests, 'interact_at', 'both');                 % 'L1'|'L2'|'both'
effect_mode         = getopt(tests, 'effect_mode', 'incremental');          % 'total'|'incremental'


[R_Xlist, R_labels] = myelin_pack(U, iR, mp);
[D_Xlist, D_labels] = myelin_pack(U, iD, mp);
[R_models, R_mlabels] = build_models_strict(R_Xlist, R_labels, modes);
[D_models, D_mlabels] = build_models_strict(D_Xlist, D_labels, modes);

% ---------------------------------------------------
% Optional debugging
if isfield(optsB,'debug') && getopt(optsB.debug,'trace_menus',false)
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


SI_all = struct(); elasticity_all = struct(); pert_all = struct();

% L1
if any(levels==1)
  % --- ROUTING ---
    base_route = struct('Xcells', {Base_route_X}, 'labels_base', {{'binary_R','caliber_R','ED'}});
    tests_R    = tests_from_models(R_models, R_mlabels);
    o = overlay_opts(optsB, struct('label_model', 'L1-Route'));
    if strcmpi(effect_mode,'incremental')                  % incremental vs total
      % incremental gain path
        if strcmpi(interaction,'none') || ismember(lower(interact_at), {'l2'})  % no ints at L1
            [SI1r, EL1r, PR1r] = bicsdual.fit_BICS_stageB2(ALLY, base_route, tests_R, pinfo, Delta_all.L1_route, o);
        else
            R_INT = apply_interactions_flagged(R_models, 'R', interaction, Base_route_X);
            INT_only = cellfun(@(main, full) full(numel(main)+1:end),R_models, R_INT, 'UniformOutput', false);
            tests_L1r = tests_from_models(INT_only, R_mlabels);
            o.reduced_models = R_models; 
            o.reduced_labels = R_mlabels; 
            o.label_model    = sprintf('L1-Route+Int{%s}', lower(interaction));
            [SI1r, EL1r, PR1r] = bicsdual.fit_BICS_stageB2(ALLY, base_route, tests_L1r, pinfo, Delta_all.L1_route, o);
        end    
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
        [SI1r, EL1r, PR1r] = bicsdual.fit_BICS_stageB2(ALLY, base_route, tests_total, pinfo, Delta_all.L1_route, o2);
        % label_out = o2.label_model;
    end
    SI_all.L1_route = SI1r; elasticity_all.L1_route = EL1r; pert_all.L1_route = PR1r;


  % --- Diffusion ---
    base_diff  = struct('Xcells', {Base_diff_X}, 'labels_base', {{'binary_D','caliber_D','ED'}});
    tests_D    = tests_from_models(D_models, D_mlabels);
    o = overlay_opts(optsB, struct('label_model', 'L1-Diff'));
    if strcmpi(effect_mode,'incremental')                  % incremental vs total
      % incremental gain path
        if strcmpi(interaction,'none') || ismember(lower(interact_at), {'l2'})  % no ints at L1
            [SI1d, EL1d, PR1d] = bicsdual.fit_BICS_stageB2(ALLY, base_diff, tests_D, pinfo, Delta_all.L1_diff, o);
        else
            D_INT = apply_interactions_flagged(D_models, 'D', interaction, Base_diff_X);
            INT_only = cellfun(@(main, full) full(numel(main)+1:end),D_models, D_INT, 'UniformOutput',false);
            tests_L1d = tests_from_models(INT_only, D_mlabels);
            o.reduced_models = D_models; 
            o.reduced_labels = D_mlabels; 
            o.label_model    = sprintf('L1-Diff+Int{%s}', lower(interaction));
            [SI1d, EL1d, PR1d] = bicsdual.fit_BICS_stageB2(ALLY, base_diff, tests_L1d, pinfo, Delta_all.L1_diff, o);
        end 
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
        [SI1d, EL1d, PR1d] = bicsdual.fit_BICS_stageB2(ALLY, base_diff, tests_total, pinfo, Delta_all.L1_diff, o2);
    end
    SI_all.L1_diff = SI1d; elasticity_all.L1_diff = EL1d; pert_all.L1_diff = PR1d;
end


% L2
if any(levels>=2)
    [RD_models, RD_labels] = cartesian_models(R_models, R_mlabels, D_models, D_mlabels);
    base_both = struct('Xcells', {Base_both_X}, 'labels_base', {{'binary_R','caliber_R','ED','binary_D','caliber_D'}});
    tests_RD  = tests_from_models(RD_models, RD_labels);
    o = overlay_opts(optsB, struct('label_model', 'L2-Both'));
    if strcmpi(effect_mode,'incremental')                  % incremental vs total
      % incremental gain path
        if strcmpi(interaction,'none') || ismember(lower(interact_at), {'l1'})  % no ints at L2
            [SI2, EL2, PR2] = bicsdual.fit_BICS_stageB2(ALLY, base_both, tests_RD, pinfo, Delta_all.L2_both, o);
        else
            RD_INT = apply_interactions_flagged(RD_models, 'RD', interaction, Base_both_X);
            INT_only = cellfun(@(main, full) full(numel(main)+1:end),RD_models, RD_INT, 'UniformOutput',false);
            tests_L2 = tests_from_models(INT_only, RD_labels);
            o.reduced_models = RD_models; 
            o.reduced_labels = RD_labels; 
            o.label_model    = sprintf('L2-Both+Int{%s}', lower(interaction));
            [SI2, EL2, PR2] = bicsdual.fit_BICS_stageB2(ALLY, base_both, tests_L2, pinfo, Delta_all.L2_both, o);
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
        [SI2, EL2, PR2] = bicsdual.fit_BICS_stageB2(ALLY, base_both, tests_total, pinfo, Delta_all.L2_both, o2);
    end
    SI_all.L2_both = SI2; elasticity_all.L2_both = EL2; pert_all.L2_both = PR2;
end

% % % % L3
% % % if any(levels==3) && doInt
% % %     [RD_models, RD_labels] = cartesian_models(R_models, R_mlabels, D_models, D_mlabels);
% % %     [INT_models, ~] = add_interactions_on_RDboth(RD_models, RD_labels, Base_both_X, {'caliber_R','caliber_D'});
% % %     tests_L3 = tests_from_models(INT_models, RD_labels); % labels same order as L2 base
% % %     base_both = struct('Xcells', {Base_both_X}, 'labels_base', {{'binary_R','caliber_R','ED','binary_D','caliber_D'}});
% % %     o = overlay_opts(optsB, struct('base_cols', numel(Base_both_X), 'reduced_models', RD_models));
% % %     [SI3, EL3, PR3] = bicsdual.fit_BICS_stageB2(ALLY, base_both, tests_L3, pinfo, Delta_all.L3_int, o);
% % %     SI_all.L3_int = SI3; elasticity_all.L3_int = EL3; pert_all.L3_int = PR3;
% % % end
% -------------------------------------------------------------------------
end
