function T = results_to_tables(stats_all, fc_labels)
%RESULTS_TO_TABLES Convert +bicsdual stats structs into tables.
levels = fieldnames(stats_all);
levels = levels(~strcmp(levels,'meta'));
global_rows = cell(0,5);
network_rows = cell(0,7);
node_rows = cell(0,6);
for li = 1:numel(levels)
    level = levels{li};
    S = stats_all.(level);
    if ~isstruct(S) || ~isfield(S,'model_labels'), continue; end
    labels = S.model_labels;
    if ischar(labels) || isstring(labels), labels = cellstr(labels); end
    nfc = numel(fc_labels);
    for m = 1:numel(labels)
        for f = 1:nfc
            fc = fc_labels{f};
            if isfield(S,'global_deltaR2')
                global_rows(end+1,:) = {level, labels{m}, fc, S.global_deltaR2(m,f), S.global_p(m,f)}; %#ok<AGROW>
            end
            if isfield(S,'ntwk_deltaR2')
                [nr, nc, ~, ~] = size(S.ntwk_deltaR2);
                if nr == numel(labels)
                    % One triangle INCLUDING within-RSN blocks. Symmetric
                    % mirrored entries are the same hypothesis, not new tests.
                    for i = 1:nc
                        for j = 1:i
                            network_rows(end+1,:) = {level, labels{m}, fc, i, j, S.ntwk_deltaR2(m,i,j,f), S.ntwk_p(m,i,j,f)}; %#ok<AGROW>
                        end
                    end
                end
            end
            if isfield(S,'node_deltaR2')
                N = size(S.node_deltaR2,2);
                for n = 1:N
                    node_rows(end+1,:) = {level, labels{m}, fc, n, S.node_deltaR2(m,n,f), S.node_p(m,n,f)}; %#ok<AGROW>
                end
            end
        end
    end
end
T = struct();
T.global = cell2table(global_rows, 'VariableNames', {'level','model','fc','deltaR2','p'});
T.network = cell2table(network_rows, 'VariableNames', {'level','model','fc','rsn_i','rsn_j','deltaR2','p'});
T.node = cell2table(node_rows, 'VariableNames', {'level','model','fc','node','deltaR2','p'});
end
