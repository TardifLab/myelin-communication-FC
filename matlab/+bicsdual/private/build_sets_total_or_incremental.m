function [RED, FULL, label] = build_sets_total_or_incremental(COL, tests)
% COL fields: base, myelin_main, int_mxcal, int_mxED  (each is a sorted unique vector of column idx)
% tests: .effect_mode {'incremental','total'}, .interactions.use_caliber, .interactions.use_ED

switch lower(tests.effect_mode)
  case 'total'
    RED  = COL.base;
    ADD  = COL.myelin_main;
    if isfield(tests,'interactions') && tests.interactions.use_caliber, ADD = union(ADD, COL.int_mxcal); end
    if isfield(tests,'interactions') && tests.interactions.use_ED,      ADD = union(ADD, COL.int_mxED);  end
    FULL = union(RED, ADD);
    label = 'Total myelin (mains+ints)';
  otherwise
    % Keep your existing incremental construction here (placeholder):
    RED  = your_current_RED;
    FULL = your_current_FULL;
    label = your_current_label;
end
end
