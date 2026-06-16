function cm = compute_communication_models(mats, cfg)
%COMPUTE_COMMUNICATION_MODELS Compute caliber, myelin, delay, and binary CMs.
%
% Output ALLX matches +bicsdual expectations:
%   {caliber, MTsat, gratio, delay, binary, ED}, where each predictor except
%   ED is a 1x6 cell array ordered as SPE, NE, SIE, PT, CMY, DE.

SCn = prep.normalize_weights(mats.caliber);
MTn = prep.normalize_weights(mats.MTsat);
GRn = prep.normalize_weights(mats.gratio);
EDn = prep.normalize_weights(mats.ED);

L_sc = prep.weight_to_length(SCn, cfg.length_transform);
L_mt = prep.weight_to_length(MTn, cfg.length_transform);
L_gr = prep.weight_to_length(GRn, cfg.length_transform);
[delay, delay_rate] = prep.make_delay(mats.length, mats.gratio, cfg.delay);

[caliber, caliber_labels, caliber_z, caliber_z_labels] = communicationModeling(SCn, L_sc, mats.ED, 'Caliber');
[mtsat,   mtsat_labels,   mtsat_z,   mtsat_z_labels]   = communicationModeling(MTn, L_mt, mats.ED, 'MTsat');
[gratio,  gratio_labels,  gratio_z,  gratio_z_labels]  = communicationModeling(GRn, L_gr, mats.ED, 'gratio');
[delaycm, delay_labels,  delay_z,   delay_z_labels]    = communicationModeling(delay_rate, delay, mats.ED, 'Delay');
[binary,  binary_labels, binary_z,  binary_z_labels]   = communicationModeling_bin(double(mats.caliber ~= 0), mats.ED, 'Binary');

cm = struct();

% Regression-ready predictors: log-transformed where appropriate and z-scored
EDz = nzzscore(EDn, 0);
cm.ALLX = {caliber_z, mtsat_z, gratio_z, delay_z, binary_z, EDz};

% Preserve raw and transformed versions for transparency.
cm.ALLX_raw = {caliber, mtsat, gratio, delaycm, binary, EDn};
cm.ALLX_z   = {caliber_z, mtsat_z, gratio_z, delay_z, binary_z, EDz};

cm.raw = struct('caliber',{caliber},'MTsat',{mtsat},'gratio',{gratio},'delay',{delaycm},'binary',{binary});
cm.z   = struct('caliber',{caliber_z},'MTsat',{mtsat_z},'gratio',{gratio_z},'delay',{delay_z},'binary',{binary_z});
cm.labels = struct('caliber',{caliber_labels},'MTsat',{mtsat_labels},'gratio',{gratio_labels},'delay',{delay_labels},'binary',{binary_labels});
cm.z_labels = struct('caliber',{caliber_z_labels},'MTsat',{mtsat_z_labels},'gratio',{gratio_z_labels},'delay',{delay_z_labels},'binary',{binary_z_labels});
cm.model_labels = cfg.communication_models;
cm.normalized_inputs = struct('caliber',SCn,'MTsat',MTn,'gratio',GRn,'ED',EDn,'delay',delay,'delay_rate',delay_rate);
end
