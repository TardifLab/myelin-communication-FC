function [Xlist, labels] = myelin_pack(U, idx_cm, predictors)
labels = {};
Xlist  = {};
for k=1:numel(predictors)
    p = lower(predictors{k});
    switch p
        case 'mtsat',   X = U.MTsat{idx_cm};  lab = 'MTsat';
        case 'gratio',  X = U.gratio{idx_cm}; lab = 'g-ratio';
        case 'delay',   X = U.delay{idx_cm};  lab = 'delay';
        otherwise, error('Unknown myelin predictor: %s', predictors{k});
    end
    Xlist{end+1} = X; labels{end+1} = sprintf('%s_%d', lab, idx_cm); %#ok<AGROW>
end
end