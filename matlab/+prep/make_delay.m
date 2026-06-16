function [delay, delay_rate] = make_delay(edge_length_mm, gratio, opts)
%MAKE_DELAY Compute simple g-ratio-based delay and inverse-delay matrices.
if nargin < 3 || isempty(opts), opts = struct; end
if ~isfield(opts,'k'), opts.k = 6.0; end
if ~isfield(opts,'axon_diameter_um'), opts.axon_diameter_um = 2.5; end
G = double(gratio);
Lmm = double(edge_length_mm);
velocity = opts.k .* (opts.axon_diameter_um ./ G);
velocity(~isfinite(velocity) | G==0) = NaN;
delay = (Lmm ./ 1000) ./ velocity;
delay(~isfinite(delay)) = 0;
delay = (delay + delay.') ./ 2;
delay(1:size(delay,1)+1:end) = 0;
delay_rate = 1 ./ delay;
delay_rate(~isfinite(delay_rate)) = 0;
delay_rate = (delay_rate + delay_rate.') ./ 2;
delay_rate(1:size(delay_rate,1)+1:end) = 0;
end
