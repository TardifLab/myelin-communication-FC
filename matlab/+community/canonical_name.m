function name = canonical_name(name)
%CANONICAL_NAME Normalize supported structural-network names.

key = regexprep(lower(char(string(name))), '[^a-z0-9]', '');
switch key
    case {'caliber','calibre','commit','commitscl','sc'}
        name = 'caliber';
    case {'mtsat','myelin','myelindensity','mc'}
        name = 'MTsat';
    case {'gratio','tractspecificgratio','mc2'}
        name = 'gratio';
    case {'delay','delays'}
        name = 'delay';
    case {'rate','delayrate','signalingrate','signallingrate'}
        name = 'rate';
    otherwise
        error('Unsupported community-detection network: %s', char(string(name)));
end
end
