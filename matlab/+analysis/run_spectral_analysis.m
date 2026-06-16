function results = run_spectral_analysis(mats, cm, pinfo, config)
%RUN_SPECTRAL_ANALYSIS Run the paper spectral-analysis workflow.
%
% Inputs
%   mats   : output of io.load_inputs
%   cm     : output of analysis.compute_communication_models
%   pinfo  : node metadata; pinfo.rsn_id is used when available
%   config : cfg.default_spectral_config(rootdir)
%
% Existing helpers reused from matlab/topology_spectra:
%   normalize_operators
%   spectral_fingerprint_windows
%   spectral_alignment_windows_robust
%   rowStandardizeWeights
%   wins_for
%   delta_topology_boxcharts (optional)

if nargin < 3 || isempty(pinfo), pinfo = struct(); end
if nargin < 4 || isempty(config)
    config = cfg.default_spectral_config();
end

spectral.validate_dependencies(config);
if ~exist(config.out_dir, 'dir'), mkdir(config.out_dir); end

inputs = spectral.get_inputs(mats, cm, config);
reference = inputs.structural.(config.reference_dataset);
N = size(reference,1);
config.kmax = min(config.kmax, N-1);
config.alignment_kmax = min(config.alignment_kmax, N-1);

if config.verbose
    fprintf('Spectral analysis: N=%d, kmax=%d\n', N, config.kmax);
end

old_visibility = get(groot, 'DefaultFigureVisible');
cleanup_visibility = onCleanup(@() set(groot, ...
    'DefaultFigureVisible', old_visibility)); %#ok<NASGU>
if config.make_plots
    set(groot, 'DefaultFigureVisible', config.figure_visible);
else
    set(groot, 'DefaultFigureVisible', 'off');
end
figures_before = findall(groot, 'Type', 'figure');

results = struct();
results.config = config;
results.inputs = rmfield(inputs, 'communication');

% Compute the caliber adjacency and normalized-Laplacian bases directly.
basis = spectral.compute_reference_basis(reference, config.kmax);
results.reference_basis = basis;

% -------------------------------------------------------------------------
% Module 1: direct pairwise spectral comparison
% -------------------------------------------------------------------------
if config.run_alignment
    if config.verbose
        fprintf('Running direct pairwise spectral alignment...\n');
    end
    alignment = struct();
    for j = 1:numel(config.target_datasets)
        name = config.target_datasets{j};
        label = config.target_labels{j};
        if config.verbose, fprintf('  caliber vs %s\n', label); end
        alignment.(name) = spectral.pairwise_alignment( ...
            reference, inputs.structural.(name), ...
            config.alignment_kmax, config.sign_align);
        if config.make_plots
            spectral.plot_pairwise_alignment(alignment.(name),label,config);
        end
    end
    results.alignment = alignment;
end

% -------------------------------------------------------------------------
% Module 2: spectral fingerprints using the existing helper
% -------------------------------------------------------------------------
if config.run_fingerprints
    if config.verbose
        fprintf('Computing communication spectral fingerprints...\n');
    end
    plot_mode = 'none';
    split_dc = false;
    if config.make_plots
        plot_mode = 'bar';
        split_dc = config.split_dc_figure;
    end

    fingerprints = spectral_fingerprint_windows( ...
        inputs.communication, inputs.model_labels, basis.Vbasis, ...
        config.windows, ...
        'Lsym', basis.Lsym, ...
        'QuantileEdges', config.quantile_edges, ...
        'StartAt', config.start_at, ...
        'RefIdx', config.reference_dataset_index, ...
        'OpForModel', config.operator_for_model, ...
        'DatasetLabels', config.dataset_labels, ...
        'Cmap', lines(numel(config.dataset_labels)), ...
        'Plot', plot_mode, ...
        'RenormCovered', config.renormalize_covered, ...
        'DeoverlapWindows', config.deoverlap_windows, ...
        'DE_UseDhalf', config.de_use_dhalf, ...
        'Dhalf', basis.Dhalf, ...
        'SplitDC', split_dc, ...
        'ShowBothDC', config.show_both_dc, ...
        'LabelMode', 'on');

    results.fingerprints = fingerprints;
    windows = fingerprints.wins;
else
    windows = spectral.default_windows(config.kmax, config.start_at);
end

% -------------------------------------------------------------------------
% Module 3: global-banded matched eigenvectors and topology metrics
% -------------------------------------------------------------------------
if config.run_matched_eigenvectors
    if config.verbose
        fprintf('Running matched-eigenvector topology analysis...\n');
    end

    target_stack = zeros(N, N, numel(config.target_datasets));
    for j = 1:numel(config.target_datasets)
        target_stack(:,:,j) = inputs.structural.(config.target_datasets{j});
    end

    opts = struct();
    opts.kmax = config.kmax;
    opts.normalize = true;
    opts.wins = windows;
    opts.doSignAlign = config.sign_align;
    opts.coords = spectral.get_optional_pinfo_field( ...
        pinfo, {'coor','coords'}, []);
    opts.SA_axis = spectral.resolve_sa_axis(config, pinfo, N);
    opts.modules = spectral.get_optional_pinfo_field( ...
        pinfo, {'rsn_id'}, []);
    if config.compute_moran
        opts.Wspatial = rowStandardizeWeights(double(reference > 0));
    else
        opts.Wspatial = [];
    end
    opts.plotMaps = false;
    opts.cmap = lines(numel(config.target_datasets));
    opts.matchMode = config.match_mode;
    opts.band = config.match_band;
    opts.LabelMode = 'on';
    opts.KPerPage = config.match_k_per_page;
    opts.dropDC = config.drop_laplacian_dc;
    opts.SArType = config.sa_corr_type;

    % The reused robust helper performs the exact global-banded matching and
    % topology calculations, but also creates large diagnostic figures.
    % Suppress those figures by default and use the lightweight plots below.
    helper_figures_before = findall(groot, 'Type', 'figure');
    visibility_before_helper = get(groot, 'DefaultFigureVisible');
    if config.make_plots && config.show_matched_helper_plots
        set(groot, 'DefaultFigureVisible', config.figure_visible);
    else
        set(groot, 'DefaultFigureVisible', 'off');
    end

    matched = spectral_alignment_windows_robust( ...
        reference, target_stack, config.target_labels, ...
        'caliber vs myelin metrics', opts);

    set(groot, 'DefaultFigureVisible', visibility_before_helper);
    helper_figures_after = findall(groot, 'Type', 'figure');
    helper_is_new = true(size(helper_figures_after));
    for hf = 1:numel(helper_figures_after)
        if ~isempty(helper_figures_before)
            helper_is_new(hf) = ~any( ...
                helper_figures_after(hf) == helper_figures_before);
        end
    end
    helper_figures = helper_figures_after(helper_is_new);
    if ~(config.make_plots && config.show_matched_helper_plots)
        close(helper_figures(ishandle(helper_figures)));
        if isfield(matched, 'figs'), matched = rmfield(matched, 'figs'); end
    end

    matched_table = spectral.matched_output_to_table(matched);
    selected_indices = spectral.load_index_selection(config, config.kmax);
    matched_table.is_selected = false(height(matched_table),1);
    is_adj = strcmp(string(matched_table.operator), 'adj');
    is_lap = strcmp(string(matched_table.operator), 'lap');
    matched_table.is_selected(is_adj) = ismember( ...
        matched_table.reference_k(is_adj), selected_indices.adj);
    matched_table.is_selected(is_lap) = ismember( ...
        matched_table.reference_k(is_lap), selected_indices.lap);

    results.matched = matched;
    results.matched_table = matched_table;
    results.selected_indices = selected_indices;

    if config.make_plots && config.plot_selected_metrics
        spectral.plot_selected_matched_metrics( ...
            matched_table(matched_table.is_selected,:), config);
    end

    if config.make_plots && config.plot_windowed_topology_boxcharts
        delta_topology_boxcharts(matched, 'adj', ...
            config.selected_metric_names, config.target_labels, ...
            lines(numel(config.target_labels)), 'on');
        delta_topology_boxcharts(matched, 'lap', ...
            config.selected_metric_names, config.target_labels, ...
            lines(numel(config.target_labels)), 'on');
    end
end

figures_after = findall(groot, 'Type', 'figure');
is_new = true(size(figures_after));
for i = 1:numel(figures_after)
    if ~isempty(figures_before)
        is_new(i) = ~any(figures_after(i) == figures_before);
    end
end
results.figure_handles = figures_after(is_new);

spectral.write_outputs(results, config);

if ~config.make_plots && ~isempty(results.figure_handles)
    close(results.figure_handles(ishandle(results.figure_handles)));
    results.figure_handles = gobjects(0);
end

if config.verbose
    fprintf('Spectral outputs written to: %s\n', config.out_dir);
end
end
