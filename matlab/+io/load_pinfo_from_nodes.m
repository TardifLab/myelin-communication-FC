function pinfo = load_pinfo_from_nodes(nodes_csv, N)
%LOAD_PINFO_FROM_NODES Build pinfo.cis and pinfo.clabels_short from nodes.csv.
%
% Required columns: node_id and either rsn_id or rsn. If node_id is absent,
% row order is used as node_id.

T = readtable(nodes_csv);
if ~ismember('node_id', T.Properties.VariableNames)
    T.node_id = (1:height(T))';
end
if height(T) ~= N
    error('nodes_csv has %d rows but matrices have %d nodes.', height(T), N);
end
if ismember('rsn_id', T.Properties.VariableNames)
    rsn_id = T.rsn_id;
else
    [~,~,rsn_id] = unique(string(T.rsn), 'stable');
end
if ismember('rsn', T.Properties.VariableNames)
    rsn = string(T.rsn);
else
    rsn = "RSN" + string(rsn_id);
end
ids = unique(rsn_id, 'stable');
pinfo = struct();
pinfo.cis = cell(1,numel(ids));
pinfo.clabels_short = cell(1,numel(ids));
for k = 1:numel(ids)
    ix = find(rsn_id == ids(k));
    pinfo.cis{k} = T.node_id(ix)';
    pinfo.clabels_short{k} = char(rsn(ix(1)));
end
pinfo.nnode = N;
pinfo.rsn_id = rsn_id;
pinfo.rsn = cellstr(rsn);
end
