function write_communication_tables(cm, mats, cfg)
%WRITE_COMMUNICATION_TABLES Save simple communication-model summary CSVs.
if ~exist(cfg.out_dir,'dir'), mkdir(cfg.out_dir); end
T = analysis.summarize_edge_correlations(cm, mats);
writetable(T, fullfile(cfg.out_dir,'communication_edge_correlations.csv'));
end
