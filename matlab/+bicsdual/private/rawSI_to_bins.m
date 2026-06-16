function Sbins = rawSI_to_bins(Smats, edges, preserveSign)
% RAWSI_TO_BINS  Map raw SI matrices to categorical integer bins.
%
% Sbins = rawSI_to_bins(Smats, edges, preserveSign)
%
% Inputs
%   Smats         : 1xK cell array of RSN×RSN raw SI matrices (z vs perm).
%   edges         : 1xN vector of nonnegative bin edges for |SI| (e.g., [2 5 10 30]).
%                   Defines N+1 bins over |SI|: [0,edges(1)), [edges(1),edges(2)), ..., [edges(N), Inf)
%   preserveSign  : (optional) logical, default true. If true, encodes sign by
%                   multiplying the bin index by sign(SI). If false, returns unsigned bins (1..N+1).
%
% Output
%   Sbins         : 1xK cell array of integer matrices with same sizes as Smats{*}.
%                   If preserveSign=true: entries are in {..., -3,-2,-1,0,+1,+2,+3,...}
%                   (0 appears only where SI==0 or NaN). If false: entries in 1..N+1.
%
% Notes
% - NaNs in SI remain NaN in Sbins.
% - Extremely large values naturally fall into the last bin ([edges(end), Inf)).
%
% Example
%   % Bin a stack of SI matrices with signed categories:
%   Sbins = rawSI_to_bins(Smats, [2 5 10 30]);       % defaults to preserveSign=true
%   % Unsigned bins (just magnitudes 1..5 for the example edges above):
%   SbinsU = rawSI_to_bins(Smats, [2 5 10 30], false);

    if nargin < 3, preserveSign = true; end
    validateattributes(Smats, {'cell'}, {'row'});
    validateattributes(edges, {'numeric'}, {'vector','nonnegative'});
    edges = sort(edges(:).');                % ensure row, ascending
    binEdges = [0, edges, inf];              % N+1 bins over |SI|

    Sbins = cell(size(Smats));
    for k = 1:numel(Smats)
        S = Smats{k};
        validateattributes(S, {'numeric'}, {'2d'});
        Sab = abs(S);

        % histcounts gives bin indices in 1..N+1 for values within [binEdges]
        % Use 'BinLimits' implicit via binEdges vector.
        [~,~,binIdx] = histcounts(Sab, binEdges);

        % Keep NaNs as NaN
        nanMask = isnan(S);
        % Optionally encode sign
        if preserveSign
            sb = sign(S) .* double(binIdx);
            % Where SI==0, sign=0 → keep 0 (distinct from bin 1)
            sb(nanMask) = NaN;
            Sbins{k} = sb;
        else
            sb = double(binIdx);
            sb(nanMask) = NaN;
            Sbins{k} = sb;
        end
    end
end
