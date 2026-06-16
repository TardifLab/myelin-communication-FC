function s = label_of(i,pinfo)
if isfield(pinfo,'clabels_short') && numel(pinfo.clabels_short)>=i && ~isempty(pinfo.clabels_short{i})
    s = pinfo.clabels_short{i};
elseif isfield(pinfo,'clabels') && numel(pinfo.clabels)>=i && ~isempty(pinfo.clabels{i})
    s = pinfo.clabels{i};
else
    s = sprintf('RSN%02d',i);
end
end
