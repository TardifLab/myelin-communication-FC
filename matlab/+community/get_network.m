function W = get_network(mats, network_name, cfg)
%GET_NETWORK Retrieve or derive the selected structural edge-weight matrix.

name = community.canonical_name(network_name);
field_name = matlab.lang.makeValidName(name);

if isfield(cfg, 'network_matrices') && isfield(cfg.network_matrices, field_name)
    W = cfg.network_matrices.(field_name);
    return
end

switch name
    case 'caliber'
        require_field(mats, 'caliber', name);
        W = mats.caliber;
    case 'MTsat'
        require_field(mats, 'MTsat', name);
        W = mats.MTsat;
    case 'gratio'
        require_field(mats, 'gratio', name);
        W = mats.gratio;
    case {'delay','rate'}
        if isfield(mats, name)
            W = mats.(name);
        else
            require_field(mats, 'length', name);
            require_field(mats, 'gratio', name);
            [delay, rate] = prep.make_delay(mats.length, mats.gratio, cfg.delay);
            if strcmp(name, 'delay'), W = delay; else, W = rate; end
        end
end

W = double(W);
if ndims(W) ~= 2 || size(W,1) ~= size(W,2)
    error('%s matrix must be square.', name);
end
end

function require_field(S, field_name, network_name)
if ~isfield(S, field_name)
    error('mats.%s is required to analyze %s.', field_name, network_name);
end
end
