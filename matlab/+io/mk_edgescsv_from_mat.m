function edges = mk_edgescsv_from_mat(config, outfile)
%MK_EDGESCSV_FROM_MAT Convert NxN matrices into a lower-triangle edge CSV.
%
% Output columns include i, j, structural predictors, one FC column per
% functional matrix, and optional rsn_i/rsn_j if node metadata are available.

if nargin < 1 || isempty(config)
    config = cfg.default_config();
end

if nargin < 2 || isempty(outfile)
    outfile = fullfile(config.out_dir, 'edges.csv');
end

[mats, pinfo] = io.load_inputs(config);

N = size(mats.caliber, 1);

% Lower triangle, excluding diagonal
[i, j] = find(tril(true(N), -1));
idx = sub2ind([N N], i, j);

edges = table();
edges.i = i;
edges.j = j;

% Structural data
edges.caliber   = mats.caliber(idx);
edges.myelin    = mats.MTsat(idx);
edges.gratio    = mats.gratio(idx);
edges.length    = mats.length(idx);
edges.euclidean = mats.ED(idx);

if isfield(mats, 'streamline_count')
    edges.streamline_count = mats.streamline_count(idx);
end

% Add all FC matrices as separate, explicitly named columns
for f = 1:numel(mats.FC)

    if isfield(config.input, 'fc_labels') && numel(config.input.fc_labels) >= f
        fc_label = matlab.lang.makeValidName(config.input.fc_labels{f});
    else
        fc_label = sprintf('%02d', f);
    end

    colname = ['FC_' fc_label];

    F = mats.FC{f};
    edges.(colname) = F(idx);
end

% Add RSN labels if node metadata are available
if isfield(pinfo, 'rsn')
    edges.rsn_i = string(pinfo.rsn(i));
    edges.rsn_j = string(pinfo.rsn(j));
end

% Remove rows with non-finite numeric values
mask = all(isfinite(edges{:, vartype('numeric')}), 2);
edges = edges(mask, :);

% Write output
outdir = fileparts(outfile);
if ~isempty(outdir) && ~exist(outdir, 'dir')
    mkdir(outdir);
end

writetable(edges, outfile);
end
