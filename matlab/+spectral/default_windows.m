function wins = default_windows(kmax,start_at)
%DEFAULT_WINDOWS Fallback windows when fingerprinting is disabled.
if nargin < 2, start_at = 2; end
wins = wins_for(kmax,'StartAt',start_at);
end
