function figs = viz_BICS_dual_allFC_summary(Delta_all, stats_all, ylbl, pinfo, opts)
% VISUALIZATION
% viz_BICS_dual without cross-FC aggregation (M=1 case only)
% - bar plot summary across groups of RSNs within FC
%       1. mean within-system (unimodal, attentional, transmodal)
%       2. mean cross-systems
%
% 2025 Mark C Nelson (MNI)
% -------------------------------------------------------------------------

if nargin<5, opts=struct; end
labelON     = getopt(opts,'labelON',true);
usek        = getopt(opts,'usek',1);
matchClim   = getopt(opts,'matchClim',false);

figpos      = [-2559 594 816 203];
fnt         = 16;
Nfc         = length(ylbl);
i_uni       = cell2mat(cellstrfind(pinfo.clabels_short,{'VIS' 'SMN'}));     % RSN groups
i_atn       = cell2mat(cellstrfind(pinfo.clabels_short,{'DAN' 'VAN'}));
i_trn       = cell2mat(cellstrfind(pinfo.clabels_short,{'CON' 'DMN'}));

labels_CMs      = stats_all.meta.label_comm_models;
labels_int      = stats_all.meta.interaction_label;

if isfield(stats_all.meta,'effect_mode') 
    labels_effect = stats_all.meta.effect_mode;
else
    labels_effect = '';
end

fields  = {'L1_route','L1_diff','L2_both','L3_int'};
nice    = {'L1-Route','L1-Diff','L2-Both','L3-Both+Int'};
f       = fields{usek};

% ====== SPECIAL CASE: all present levels have M==1 ======
figs = struct();
mean_uni = nan(1,Nfc);
mean_atn = nan(1,Nfc);
mean_trn = nan(1,Nfc);

% ---- Gather data
for k=1:Nfc
    
    St = stats_all.(f);
    blk_ij = St.block_map; Nntwk = St.Nntwk;

    % FC
    Dt  = Delta_all.(f)(:,1,k);    % B x 1 x Nfc -> B x 1

    % Map to Nntwk x Nntwk matrices
      Dm = blocks_to_mat(Dt(:,1), blk_ij, Nntwk, true);

    % Masks for within-system
      mask_uni_vec  = all(ismember(blk_ij, i_uni),2);
      mask_atn_vec  = all(ismember(blk_ij, i_atn),2);
      mask_trn_vec  = all(ismember(blk_ij, i_trn),2);
      mask_uni_mat  = logical(blocks_to_mat(mask_uni_vec, blk_ij, Nntwk, true));
      mask_atn_mat  = logical(blocks_to_mat(mask_atn_vec, blk_ij, Nntwk, true));
      mask_trn_mat  = logical(blocks_to_mat(mask_trn_vec, blk_ij, Nntwk, true));
      mask_lt       = logical(tril(ones(Nntwk)));

    % --- Masks for cross-system
    % Helper inline: true if edge connects set A <-> set B (either direction)
      is_cross = @(A,B) ( ismember(blk_ij(:,1),A) & ismember(blk_ij(:,2),B) ) | ...
                        ( ismember(blk_ij(:,1),B) & ismember(blk_ij(:,2),A) );
    % Cross-system masks
      mask_uni_atn_vec  = is_cross(i_uni, i_atn);
      mask_atn_trn_vec  = is_cross(i_atn, i_trn);
      mask_uni_trn_vec  = is_cross(i_uni, i_trn);
      mask_uni_atn_mat  = logical(blocks_to_mat(mask_uni_atn_vec, blk_ij, Nntwk, true));
      mask_atn_trn_mat  = logical(blocks_to_mat(mask_atn_trn_vec, blk_ij, Nntwk, true));
      mask_uni_trn_mat  = logical(blocks_to_mat(mask_uni_trn_vec, blk_ij, Nntwk, true));
      

    % Summarize
      mean_uni(k)       = mean(Dm(mask_uni_mat & mask_lt));
      mean_atn(k)       = mean(Dm(mask_atn_mat & mask_lt));
      mean_trn(k)       = mean(Dm(mask_trn_mat & mask_lt));
      mean_uni_atn(k)   = mean(Dm(mask_uni_atn_mat & mask_lt));
      mean_atn_trn(k)   = mean(Dm(mask_atn_trn_mat & mask_lt));
      mean_uni_trn(k)   = mean(Dm(mask_uni_trn_mat & mask_lt));
end

DPLT    = [mean_uni; mean_atn; mean_trn];
cmaptmp = gray(3+1);


% ---- Figure A: bar plot within-system
figs.barplot = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect ', ' nice{usek} ', within-system'],figpos); hold on
  B = bar(DPLT', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  for k=1:size(DPLT,1), set(B(k),'FaceColor',cmaptmp(k,:)); end  
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fnt+2,'Helvetica'); 
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('\\mu %s','\DeltaR^2')); 
  set(gca,'TickLength',[0 0]); % legend({'unimodal' 'attentional' 'transmodal'})


DPLT    = [mean_uni_atn; mean_atn_trn; mean_uni_trn];

% ---- Figure B: bar plot cross-system
figs.barplot2 = myfig([labels_CMs ', ' labels_int ', effect-' labels_effect ', ' nice{usek} ', cross-system'],figpos); hold on
  B = bar(DPLT', 'grouped'); grid on; set(gca,'XGrid','off'); box on;
  for k=1:size(DPLT,1), set(B(k),'FaceColor',cmaptmp(k,:)); end  
  xticks(1:Nfc); xticklabels(ylbl); xtickangle(30); font(fnt+2,'Helvetica'); 
  for k=1:Nfc-1; line([k+.5 k+.5],ylim,'Color','k','LineStyle','-','LineWidth',0.5); end   
  ylabel(sprintf('\\mu %s','\DeltaR^2')); 
  set(gca,'TickLength',[0 0]); % legend({'uni-atn' 'atn-trn' 'uni-trn'})




% -------------------------------------------------------------------------
end
