function [RED, FULL, effect_label] = build_sets_effect_mode(effect_mode, interaction, COL, RED_inc, FULL_inc, label_inc)
% BUILD_SETS_EFFECT_MODE
% Switch between incremental (existing behavior) vs total myelin contribution.
%
% Inputs
%   effect_mode : char  'incremental' | 'total'
%   interaction : char  'none' | 'caliber' | 'ed' | 'both'
%   COL         : struct with integer index vectors for current test *scope*
%                 .base          -> base columns (e.g., binary_cmi, caliber_cmi, ED)
%                 .myelin_main   -> myelin main-effect columns in-scope
%                 .int_caliber   -> myelin×caliber interaction columns in-scope (or [])
%                 .int_ed        -> myelin×ED      interaction columns in-scope (or [])
%                 (Important: COL must already respect 'interact_at' and Level (L1/L2))
%   RED_inc     : integer vector (your existing incremental RED)
%   FULL_inc    : integer vector (your existing incremental FULL)
%   label_inc   : char (your existing label for incremental tests; optional)
%
% Outputs
%   RED, FULL   : integer vectors to pass to the model fit
%   effect_label: char, human-readable label for reporting/plots
%
% Notes
% - Does not alter existing logic for how you currently build incremental RED/FULL.
% - For 'total', it compares BASE vs BASE+MYELIN(mains)+requested interactions.
%
% Mark C Nelson / dual BICS extension (2025-11-06)

if nargin < 6 || isempty(label_inc), label_inc = 'incremental'; end

% Basic sanity (avoid accidental overlaps that can inflate df)
assert(isvector(COL.base)        && isnumeric(COL.base),        'COL.base must be a numeric vector');
assert(isvector(COL.myelin_main) && isnumeric(COL.myelin_main), 'COL.myelin_main must be a numeric vector');
if ~isfield(COL,'int_caliber'), COL.int_caliber = []; end
if ~isfield(COL,'int_ed'),      COL.int_ed      = []; end

% Ensure uniqueness/sorting defensively
COL.base        = unique(COL.base(:))';
COL.myelin_main = unique(COL.myelin_main(:))';
COL.int_caliber = unique(COL.int_caliber(:))';
COL.int_ed      = unique(COL.int_ed(:))';

switch lower(effect_mode)
  case 'incremental'
    % Preserve current behavior exactly
    RED          = unique(RED_inc(:))';
    FULL         = unique(FULL_inc(:))';
    effect_label = label_inc;

  case 'total'
    % Total myelin contribution: FULL = BASE ∪ MAIN ∪ (requested interactions), RED = BASE
    RED = COL.base;

    ADD = COL.myelin_main;
    switch lower(interaction)
      case 'none'
        % add no interactions
      case 'caliber'
        ADD = union(ADD, COL.int_caliber);
      case 'ed'
        ADD = union(ADD, COL.int_ed);
      case 'both'
        ADD = union(union(ADD, COL.int_caliber), COL.int_ed);
      otherwise
        error('Unknown interaction option: %s', interaction);
    end

    FULL = union(RED, ADD);

    % Label for reporting
    switch lower(interaction)
      case 'none',  i_tag = 'mains only';
      case 'caliber', i_tag = 'mains + myelin×caliber';
      case 'ed',      i_tag = 'mains + myelin×ED';
      case 'both',    i_tag = 'mains + both interactions';
    end
    effect_label = sprintf('Total myelin (%s)', i_tag);

  otherwise
    error('Unknown effect_mode: %s (use ''incremental'' or ''total'')', effect_mode);
end

% Final safety: avoid empty FULL/RED or degenerate cases
assert(~isempty(RED)  && all(RED  > 0), 'RED set is empty or invalid.');
assert(~isempty(FULL) && all(FULL > 0), 'FULL set is empty or invalid.');
assert(all(ismember(RED,  FULL)) || isempty(intersect(FULL, setdiff(1:max([FULL,RED]), union(FULL,RED)))), ...
  'Unexpected: RED not subset of FULL in total mode or column sets malformed.');

end
