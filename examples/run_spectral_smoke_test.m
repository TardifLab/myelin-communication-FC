%% Synthetic smoke test for the spectral add-on
rootdir = fileparts(fileparts(mfilename('fullpath')));
cd(fullfile(rootdir,'matlab'));
setup_paths(rootdir);
rng(11);

N = 40;
mask = triu(rand(N)<0.22,1);
mask = mask + mask.';
makeW = @() local_weighted(mask,N);

mats = struct();
mats.caliber = makeW();
mats.MTsat = makeW();
mats.gratio = makeW();
mats.length = 20 + 80*rand(N);
mats.length = (mats.length + mats.length')./2;
mats.length(~mask) = 0;
mats.length(1:N+1:end) = 0;
mats.ED = rand(N);
mats.ED = (mats.ED + mats.ED')./2;
mats.ED(1:N+1:end) = 0;

cm = struct();
cm.model_labels = {'SPE','NE','SIE','PT','CMY','DE'};
for nm = {'caliber','MTsat','gratio','delay'}
    C = cell(1,6);
    for m = 1:6
        X = rand(N);
        X = (X + X')./2;
        X(1:N+1:end) = 0;
        C{m} = X;
    end
    cm.raw.(nm{1}) = C;
end
cm.normalized_inputs.caliber = mats.caliber;
cm.normalized_inputs.MTsat = mats.MTsat;
cm.normalized_inputs.gratio = mats.gratio;
cm.normalized_inputs.delay_rate = makeW();

pinfo.rsn_id = repelem((1:4)',N/4);

config = cfg.default_spectral_config(rootdir);
config.out_dir = fullfile(rootdir,'results','spectral_smoke_test');
config.kmax = N-1;
config.alignment_kmax = 20;
config.make_plots = false;
config.matched_indices_mode = 'custom';
config.custom_indices.adj = 2:10;
config.custom_indices.lap = 2:10;
results = analysis.run_spectral_analysis(mats,cm,pinfo,config); %#ok<NASGU>
fprintf('Spectral smoke test completed.\n');

function W = local_weighted(mask,N)
W = rand(N).*mask;
W = (W+W')./2;
W(1:N+1:end) = 0;
end
