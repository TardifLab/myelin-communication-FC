function results = run_community_detection(mats, config)
%RUN_COMMUNITY_DETECTION Load or rerun structural-network community analysis.
%
% Precomputed use:
%   cfg = cfg.default_community_config();
%   results = analysis.run_community_detection([],cfg);
%
% Rerun use:
%   [mats,~] = io.load_inputs(main_cfg);
%   cfg.mode = 'rerun';
%   cfg.selection_mode = 'none';
%   results = analysis.run_community_detection(mats,cfg);

if nargin < 2 || isempty(config)
    config = cfg.default_community_config();
end
if nargin < 1, mats = []; end
community.validate_config(config);

mode = lower(char(string(config.mode)));
if ischar(config.networks) || isstring(config.networks)
    networks = cellstr(config.networks);
else
    networks = config.networks;
end
if ~exist(config.out_dir,'dir'), mkdir(config.out_dir); end
results = struct();

for k = 1:numel(networks)
    name = community.canonical_name(networks{k});
    fprintf('\nCommunity detection: %s (%s mode)\n',name,mode);

    switch mode
        case 'precomputed'
            R = community.load_precomputed(name,config.precomputed_dir);
        case 'rerun'
            if isempty(mats)
                error('Rerun mode requires the loaded mats structure.');
            end
            if strcmpi(config.selection_mode,'paper')
                error(['selection_mode = ''paper'' is available only in precomputed mode. ', ...
                    'Use ''none'', ''score'', or ''gamma'' when rerunning.']);
            end
            W = community.get_network(mats,name,config);
            R = community.run_two_pass(W,name,config);
        otherwise
            error('cfg.mode must be ''precomputed'' or ''rerun''.');
    end

    selected = community.select_partition(R,config);
    paths = community.write_outputs(R,selected,config);

    if config.plot_diagnostics
        H = community.plot_stability(R,config,selected);
        if config.save_figures
            figdir = fullfile(config.out_dir,lower(name));
            paths.diagnostics_figure = fullfile(figdir,'stability_diagnostics.png');
            exportgraphics(H,paths.diagnostics_figure,'Resolution',200);
        end
    end

    field_name = matlab.lang.makeValidName(name);
    results.(field_name) = struct('sweep',R,'selected',selected,'paths',paths);
end
end
