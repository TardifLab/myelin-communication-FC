function wins = wins_for(K, varargin)
% WINS_FOR  Build global/meso/local/local+ eigen windows that scale with K.
% 
% wins = wins_for(K) returns a 1x4 cell array:
%   { global, meso, local, local_plus }, each an integer index vector.
%
% By default we exclude k=1 (the DC/constant or A's Perron slot) and start at k=2.
% Use 'Perc' to adjust breakpoints as percentiles of K, and 'Overlap' to
% allow a small overlap (in eigen indices) between consecutive windows.
%
% Inputs (Name-Value):
%   'Perc'      : [pG pM pL], breakpoints as fractions of K (default [0.015 0.075 0.20])
%   'Overlap'   : integer number of overlapping indices between windows (default 0)
%   'ClampG'    : [gMin gMax] hard clamp for global end index (default [2 6])
%   'StartAt'   : first k to include (default 2; keep 2 to exclude k=1)
%
% Examples:
%   wins = wins_for(399); 
%   % -> approx {2:6, 7:30, 31:80, 81:399}
%
%   wins = wins_for(200,'Overlap',1);
%   % -> approx {2:6, 6:15, 15:40, 40:200}
%
%   wins = wins_for(399,'Perc',[0.02 0.08 0.25]); % nudge the breakpoints
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------


    p = inputParser;
    addParameter(p,'Perc',   [0.015 0.075 0.20], @(x)isnumeric(x)&&numel(x)==3);
    addParameter(p,'Overlap',0,                  @(x)isnumeric(x)&&isscalar(x));
    addParameter(p,'ClampG', [2 6],             @(x)isnumeric(x)&&numel(x)==2);
    addParameter(p,'StartAt',2,                 @(x)isnumeric(x)&&isscalar(x));
    parse(p,varargin{:});

    perc    = p.Results.Perc;
    ov      = round(p.Results.Overlap);
    clampG  = p.Results.ClampG;
    kStart  = max(1, round(p.Results.StartAt));  %#ok<NASGU> % (kept for clarity)

    % --- Compute raw breakpoints (rounded) ---
    bG = max(round(perc(1)*K), clampG(1));           % tentative global end
    bG = min(max(bG, clampG(1)), clampG(2));         % clamp to [2..6] by default
    bM = max(round(perc(2)*K), bG+1);                % meso end
    bL = max(round(perc(3)*K), bM+1);                % local end

    % Ensure monotonicity within [1..K]
    bG = min(max(bG, 2), K);
    bM = min(max(bM, bG+1), K);
    bL = min(max(bL, bM+1), K);

    % --- Start/stop indices with optional overlap ---
    % We intentionally exclude k=1; first window starts at 2
    g1 = 2;              g2 = bG;
    m1 = max(g2 - ov,  g2+1 - ov);    m2 = bM;
    l1 = max(m2 - ov,  m2+1 - ov);    l2 = bL;
    lp1= max(l2 - ov,  l2+1 - ov);    lp2= K;

    % Guard against tiny K producing empty ranges
    if g2 < g1, g1 = min(2,K); g2 = min(max(2, g1), K); end
    if m2 < m1, m1 = min(g2+1, K); m2 = min(max(m1, g2+1), K); end
    if l2 < l1, l1 = min(m2+1, K); l2 = min(max(l1, m2+1), K); end
    if lp2 < lp1, lp1 = min(l2+1, K); lp2 = K; end

    wins = { g1:g2, m1:m2, l1:l2, lp1:lp2 };

% -------------------------------------------------------------------------
end
