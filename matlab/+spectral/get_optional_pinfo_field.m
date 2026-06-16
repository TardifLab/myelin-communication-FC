function value = get_optional_pinfo_field(pinfo, candidates, default_value)
%GET_OPTIONAL_PINFO_FIELD Return the first available pinfo field.
value = default_value;
if isempty(pinfo) || ~isstruct(pinfo), return; end
for i = 1:numel(candidates)
    if isfield(pinfo,candidates{i}) && ~isempty(pinfo.(candidates{i}))
        value = pinfo.(candidates{i});
        return
    end
end
end
