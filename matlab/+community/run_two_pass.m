function R = run_two_pass(W, network_name, cfg)
%RUN_TWO_PASS Run coarse-to-fine Louvain community detection.

community.require_dependencies();
name = community.canonical_name(network_name);
[Dxfm, prep_meta] = community.preprocess_network(W, cfg);
N = size(Dxfm,1);

% ---------------- First pass: coarse gamma range ----------------
gamma_coarse = logspace(log10(cfg.coarse.gamma_bounds(1)), ...
    log10(cfg.coarse.gamma_bounds(2)), cfg.coarse.n_gamma).';
coarse_partitions = zeros(N, cfg.coarse.n_gamma, 'uint16');
coarse_q = nan(cfg.coarse.n_gamma,1);
coarse_n = zeros(cfg.coarse.n_gamma,1);

for g = 1:cfg.coarse.n_gamma
    rng(cfg.random_seed + g - 1, 'twister');
    [ci,q] = community_louvain(Dxfm, gamma_coarse(g), [], []);
    coarse_partitions(:,g) = uint16(ci);
    coarse_q(g) = q;
    coarse_n(g) = numel(unique(ci));
end

Yu = find(coarse_n > cfg.coarse.lower_count_threshold, 1, 'first');
Yv = find(coarse_n > (N*cfg.coarse.upper_count_fraction + 1), 1, 'first');
if isempty(Yu), Yu = 1; end
if isempty(Yv), Yv = numel(gamma_coarse); end
if Yu >= Yv
    Yu = max(1, floor(0.25*numel(gamma_coarse)));
    Yv = min(numel(gamma_coarse), ceil(0.90*numel(gamma_coarse)));
end

pad = cfg.fine.padding;
gamma_min = gamma_coarse(max(1,Yu-pad));
gamma_max = gamma_coarse(min(numel(gamma_coarse),Yv+pad));
gamma = logspace(log10(gamma_min), log10(gamma_max), cfg.fine.n_gamma).';

% ---------------- Second pass: repeated partitions and consensus --------
consensus_partitions = zeros(N, cfg.fine.n_gamma, 'uint16');
modularity = nan(cfg.fine.n_gamma, cfg.fine.n_repetitions);
zrand_mean = nan(cfg.fine.n_gamma,1);
zrand_variance = nan(cfg.fine.n_gamma,1);

for g = 1:cfg.fine.n_gamma
    fprintf('%s: gamma %d/%d\n', name, g, cfg.fine.n_gamma);
    ci_reps = zeros(N, cfg.fine.n_repetitions, 'uint16');
    q_reps = nan(1, cfg.fine.n_repetitions);
    seeds = cfg.random_seed + 10000*g + (1:cfg.fine.n_repetitions);

    if cfg.use_parallel
        parfor r = 1:cfg.fine.n_repetitions
            rng(seeds(r),'twister');
            [ci,q] = community_louvain(Dxfm, gamma(g), [], []);
            ci_reps(:,r) = uint16(ci);
            q_reps(r) = q;
        end
    else
        for r = 1:cfg.fine.n_repetitions
            rng(seeds(r),'twister');
            [ci,q] = community_louvain(Dxfm, gamma(g), [], []);
            ci_reps(:,r) = uint16(ci);
            q_reps(r) = q;
        end
    end
    modularity(g,:) = q_reps;

    ag = agreement(double(ci_reps)) ./ cfg.fine.n_repetitions;
    ci_null = zeros(size(ci_reps),'uint16');
    for r = 1:cfg.fine.n_repetitions
        rng(cfg.random_seed + 200000 + 10000*g + r, 'twister');
        ci_null(:,r) = ci_reps(randperm(N),r);
    end
    ag_null = agreement(double(ci_null)) ./ cfg.fine.n_repetitions;
    tau = mean(ag_null(:));
    rng(cfg.random_seed + 300000 + g, 'twister');
    consensus_partitions(:,g) = uint16(consensus_und(ag, tau, cfg.fine.consensus_repetitions));

    zvals = nan(cfg.fine.n_repetitions*(cfg.fine.n_repetitions-1)/2,1);
    c = 0;
    for i = 1:(cfg.fine.n_repetitions-1)
        for j = (i+1):cfg.fine.n_repetitions
            c = c + 1;
            zvals(c) = community.zrand(ci_reps(:,i), ci_reps(:,j));
        end
    end
    zrand_mean(g) = mean(zvals,'omitnan');
    zrand_variance(g) = var(zvals,0,'omitnan');
end

R = struct();
R.schema_version = '1.0';
R.network_name = name;
R.source_label = name;
R.n_nodes = N;
R.gamma = gamma;
R.consensus_partitions = consensus_partitions;
R.n_communities = arrayfun(@(k) numel(unique(consensus_partitions(:,k))), ...
    (1:size(consensus_partitions,2))).';
R.zrand_mean = zrand_mean;
R.zrand_variance = zrand_variance;
R.zrand_valid = isfinite(zrand_mean) & isfinite(zrand_variance);
R.zrand_mean_raw = zrand_mean;
R.modularity = modularity;
R.mean_modularity = mean(modularity,2,'omitnan');
R.input_matrix = double(W);
R.transformed_matrix = Dxfm;
R.paper = struct();
R.mode = 'rerun';
R.preprocess = prep_meta;
R.coarse = struct('gamma',gamma_coarse,'partitions',coarse_partitions, ...
    'modularity',coarse_q,'n_communities',coarse_n, ...
    'lower_index',Yu,'upper_index',Yv, ...
    'fine_gamma_bounds',[gamma_min gamma_max]);
end
