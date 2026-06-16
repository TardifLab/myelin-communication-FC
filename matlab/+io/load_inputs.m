function [mats, pinfo] = load_inputs(cfg)
%LOAD_INPUTS Load matrices and parcellation/network metadata from cfg.

mats = struct();
mats.caliber  = io.load_matrix_var(cfg.input.caliber_mat);
mats.MTsat    = io.load_matrix_var(cfg.input.mtsat_mat);
mats.gratio   = io.load_matrix_var(cfg.input.gratio_mat);
mats.length   = io.load_matrix_var(cfg.input.length_mat);
mats.ED       = io.load_matrix_var(cfg.input.euclidean_mat);

nfc = numel(cfg.input.fc_mats);
mats.FC = cell(1,nfc);
for f = 1:nfc
    mats.FC{f} = io.load_matrix_var(cfg.input.fc_mats{f});
end
mats.fc_labels = cfg.input.fc_labels;

N = size(mats.caliber,1);
if isfield(cfg.input,'nodes_csv') && exist(cfg.input.nodes_csv,'file')
    pinfo = io.load_pinfo_from_nodes(cfg.input.nodes_csv, N);
else
    pinfo = io.default_pinfo(N);
end
end
