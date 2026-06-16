function L = ensure_labels(s, field, fallback)
L = fallback;
if isstruct(s) && isfield(s, field) && ~isempty(s.(field))
    L = s.(field);
end
end
