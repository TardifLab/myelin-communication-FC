function files = write_fdr_tables(F,outdir)
%WRITE_FDR_TABLES Export each unique model once, with raw and adjusted values.
% deltaR2 and p remain raw; deltaR2_fdr is a NaN-masked display copy.
% Existing exports are protected. Fits/spin tests are not run here.
if ~exist(outdir,'dir'), mkdir(outdir); end
refs = [1 1;3 1;1 2;2 2;1 3;2 3;3 3;4 3];
levels = {'L1_route','L1_diff','L2_both'};
scales = {'global','network','node'};
files = cell(numel(F.meta.conditions),numel(scales));
% Check all destinations before writing any table.
for c = 1:numel(F.meta.conditions)
    for s = 1:numel(scales)
        files{c,s} = fullfile(outdir,sprintf('BICS_FDR_%s_%s.csv', ...
            scales{s},F.meta.conditions{c}));
        assert(~isfile(files{c,s}),'BICSFDR:OutputExists', ...
            'Output exists: %s. Choose a new output directory.',files{c,s});
    end
end
for c = 1:numel(F.meta.conditions)
    cond = F.meta.conditions{c};
    combined = cell(1,numel(scales));
    for m = 1:size(refs,1)
        pair = refs(m,1); lev = levels{refs(m,2)};
        A = F.pairs(pair).conditions.(cond);
        one = struct(); one.(lev) = A.stats_all.(lev);
        tables = analysis.results_to_tables(one,F.meta.FCLabels);
        st = one.(lev);
        for s = 1:numel(scales)
            scale = scales{s}; T = tables.(scale);
            if isempty(T), continue; end
            [found,fc_index] = ismember(string(T.fc),string(F.meta.FCLabels));
            assert(all(found),'BICSFDR:Labels','FDR table FC labels do not match.');
            q = nan(height(T),1);
            for row = 1:height(T)
                f = fc_index(row);
                switch scale
                    case 'global'
                        q(row) = st.global_q(1,f);
                    case 'network'
                        q(row) = st.ntwk_q(1,T.rsn_i(row),T.rsn_j(row),f);
                    case 'node'
                        q(row) = st.node_q(1,T.node(row),f);
                end
            end
            T.q = q;
            T.reject_fdr = isfinite(q) & q<F.meta.alpha;
            T.deltaR2_fdr = T.deltaR2;
            T.deltaR2_fdr(~T.reject_fdr) = NaN;
            T.communication_model = repmat(string(F.meta.model_order{m}),height(T),1);
            T.fdr_family = repmat(string(st.fdr_family),height(T),1);
            T.condition = repmat(string(cond),height(T),1);
            T.fdr_method = repmat("BH",height(T),1);
            T.fdr_alpha = repmat(F.meta.alpha,height(T),1);
            if isempty(combined{s}), combined{s} = T;
            else, combined{s} = [combined{s};T]; end %#ok<AGROW>
        end
    end
    for s = 1:numel(scales)
        assert(~isempty(combined{s}),'BICSFDR:EmptyTable','Missing %s data.',scales{s});
        tmp = [tempname(outdir) '.csv'];
        cleanup = onCleanup(@()delete_temp(tmp));
        writetable(combined{s},tmp);
        [ok,msg] = movefile(tmp,files{c,s});
        assert(ok,'BICSFDR:Save','%s',msg);
        clear cleanup
    end
end
end

function delete_temp(path)
if isfile(path), delete(path); end
end
