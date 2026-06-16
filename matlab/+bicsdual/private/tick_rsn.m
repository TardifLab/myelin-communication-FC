function tick_rsn(pinfo,Nntwk)
xticks(1:Nntwk); yticks(1:Nntwk);
if isfield(pinfo,'clabels_short') && numel(pinfo.clabels_short)==Nntwk
    xticklabels(pinfo.clabels_short); yticklabels(pinfo.clabels_short);
    xtickangle(40); ytickangle(0);
end
end