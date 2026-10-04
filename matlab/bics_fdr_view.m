function [Delta_all, stats_all, Q_all, H_all] = bics_fdr_view(F, pairIndex, condition, mode)
%BICS_FDR_VIEW Return independent display copies; never modifies F.
% [Delta_all,stats_all,Q_all,H_all] = bics_fdr_view(F,1,'total','fdr');
% Modes: 'raw' (default), 'fdr', 'uncorrected'.
% 'fdr' masks estimates whose adjusted p is missing or >= F.meta.alpha.
% 'uncorrected' uses raw p at the same alpha, for comparison only.
% Masking uses NaN, NOT zero. All p/q fields retain their original values.
% meta.selected_condition and meta.display_mode describe this view. The
% original interaction_label is preserved: int-both can also describe total.
% Use RAW stats_all for S-A correlations and descriptive block averages.
% NaN display colors depend on the plotting function. bics_fdr_preview
% explicitly renders missing/masked blocks gray using AlphaData.
% Do not pass fdr-masked stats into sa_corr_nodal for inference: it would
% calculate correlations on the selected nodes. Extract node_deltaR2 only
% for surface display, and use raw stats for subsequent statistical work.
if nargin<4, mode='raw'; end
assert(isscalar(pairIndex)&&pairIndex>=1&&pairIndex<=numel(F.pairs)&&pairIndex==fix(pairIndex), ...
    'BICSFDR:PairIndex','Invalid pair index.');
condition=char(condition); mode=lower(char(mode));
assert(ismember(mode,{'raw','fdr','uncorrected'}),'BICSFDR:ViewMode','Use raw, fdr, or uncorrected.');
assert(isfield(F.pairs(pairIndex).conditions,condition),'BICSFDR:ViewCondition','Condition not found.');
A=F.pairs(pairIndex).conditions.(condition);
Delta_all=A.Delta_all; stats_all=A.stats_all; Q_all=A.Q_all; H_all=A.H_all;
if ~isfield(stats_all,'meta'), stats_all.meta=struct; end
stats_all.meta.selected_condition=condition;
stats_all.meta.display_mode=mode;
if strcmp(mode,'raw'), return; end
levels={'L1_route','L1_diff','L2_both'};
for k=1:3
    lev=levels{k}; st=stats_all.(lev);
    if strcmp(mode,'fdr')
        h=A.H_all.(lev); hg=st.global_h; hn=st.node_h; ht=st.ntwk_h;
    else
        p=A.P_all.(lev); h=isfinite(p)&p<F.meta.alpha;
        hg=isfinite(st.global_p)&st.global_p<F.meta.alpha;
        hn=isfinite(st.node_p)&st.node_p<F.meta.alpha;
        ht=isfinite(st.ntwk_p)&st.ntwk_p<F.meta.alpha;
    end
    d=Delta_all.(lev); d(~h)=NaN; Delta_all.(lev)=d;
    st.global_deltaR2(~hg)=NaN;
    st.node_deltaR2(~hn)=NaN;
    st.ntwk_deltaR2(~ht)=NaN;
    stats_all.(lev)=st;
end
end
