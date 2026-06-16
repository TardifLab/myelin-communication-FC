function T = summarize_edge_correlations(cm, mats)
%SUMMARIZE_EDGE_CORRELATIONS Spearman correlations among communication matrices.
%
% Correlations are summarized across:
%   all = all lower-triangle node pairs
%   con = node pairs with an empirical structural connection
%   dis = node pairs without an empirical structural connection
%
% Connected/disconnected masks are derived from the empirical binary
% structural network, here defined from nonzero caliber/MTsat/gratio edges.

sets = {'caliber','MTsat','gratio','delay','binary'};
model_labels = cm.model_labels;

N = size(mats.caliber,1);
lt_mask = tril(true(N), -1);

% Empirical connection mask.
% Use caliber/MTsat/gratio jointly so the mask is not dependent on one
% particular microstructural edge definition.
empirical_con = (mats.caliber ~= 0) | (mats.MTsat ~= 0) | (mats.gratio ~= 0);
empirical_con = empirical_con & lt_mask;
empirical_dis = ~empirical_con & lt_mask;

edge_groups = struct();
edge_groups.all = lt_mask;
edge_groups.con = empirical_con;
edge_groups.dis = empirical_dis;

group_names = {'all','con','dis'};

rows = {};

for m = 1:numel(model_labels)
    for a = 1:numel(sets)
        for b = a+1:numel(sets)

            Xa = cm.raw.(sets{a}){m};
            Xb = cm.raw.(sets{b}){m};

            for g = 1:numel(group_names)

                group = group_names{g};
                mask = edge_groups.(group);

                x = Xa(mask);
                y = Xb(mask);

                ok = isfinite(x) & isfinite(y);

                n_edges = sum(ok);

                if n_edges < 3 || numel(unique(x(ok))) < 2 || numel(unique(y(ok))) < 2
                    rho = NaN;
                else
                    rho = corr(x(ok), y(ok), ...
                        'type', 'Spearman', ...
                        'rows', 'complete');
                end

                rows(end+1,:) = { ...
                    model_labels{m}, ...
                    sets{a}, ...
                    sets{b}, ...
                    group, ...
                    rho, ...
                    n_edges}; %#ok<AGROW>
            end
        end
    end
end

T = cell2table(rows, ...
    'VariableNames', {'communication_model','predictor_a','predictor_b','edge_group','spearman_rho','n_edges'});
end
