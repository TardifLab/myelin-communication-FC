function Y_out = transform_fc(Y_in, ylbl)
% TRANSFORM_FC
% Apply appropriate variance-stabilizing transforms to each FC matrix
% based on the dataset label. Keeps matrices symmetric with zero diagonals.
%
% Heuristics:
%   - 'bold', 'corr', 'aec' -> Fisher r-to-z (atanh), clipped at ±(1-eps)
%   - 'coh' (magnitude-squared coherence), 'plv', 'wpli', 'dwpli' -> logit on [0,1]
%   - 'icoh' (imaginary coherence, [-1,1]) -> atanh on [-1,1], but interpret cautiously
%   - otherwise -> leave as-is
%
%
% 2025 Mark C Nelson, McConnell Brain Imaging Centre, MNI, McGill
%--------------------------------------------------------------------------


    assert(numel(Y_in)==numel(ylbl), 'Y_in and ylbl must match in length.');
    Y_out = Y_in;
    for d = 1:numel(Y_in)
        A = Y_in{d};
        lbl = lower(string(ylbl{d}));

        A = (A + A.')/2;                  % enforce symmetry
        A(1:size(A,1)+1:end) = 0;         % clean diagonal

        if contains(lbl,'bold') || contains(lbl,'corr') || contains(lbl,'aec')
            % Correlation-like metrics (Pearson of envelopes, etc.)
            A = atanh( clip(A, -0.999999, 0.999999) );

        elseif contains(lbl,'icoh') || contains(lbl,'imag')
            % Imaginary coherence can be negative; bounded (-1,1)
            % atanh is defined, but sampling distribution differs from Pearson’s r.
            % Use with caution; alternatives: rank-normalize or null z-scoring.
            A = atanh( clip(A, -0.999999, 0.999999) );

        elseif contains(lbl,'coh') || contains(lbl,'plv') || contains(lbl,'wpli')
            % [0,1]-bounded measures; use logit
            A = logit( clip(A, 1e-6, 1-1e-6) );

        else
            % Unknown metric: leave as-is (or choose a conservative rank-normalization)
            A = A;
        end

        % Re-set diagonal to 0 after transform (keeps block means clean)
        A(1:size(A,1)+1:end) = 0;
        Y_out{d} = A;
    end
end

function X = clip(X, lo, hi)
    X = min(max(X, lo), hi);
end

function Z = logit(P)
    Z = log(P ./ (1 - P));
end
