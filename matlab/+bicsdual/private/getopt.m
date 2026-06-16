function v = getopt(s, name, default)
if ~isstruct(s), v = default; return; end
if isfield(s,name), v=s.(name); else, v=default; end
end

% % % function v = getopt(s, name, default)
% % % % Robust getopt:
% % % % - Works with scalar structs and struct arrays (captures comma-separated lists)
% % % % - Preserves cells; converts string arrays to cellstr; leaves numerics alone
% % % 
% % % v = default;
% % % if ~isstruct(s) || ~isfield(s,name), return; end
% % % 
% % % val = s.(name);
% % % 
% % % % If 's' is a struct array, s.(name) yields a comma-separated list.
% % % % Wrap it to capture all elements into a cell array.
% % % if numel(s) > 1
% % %     val = {s.(name)};  % collect comma-separated list
% % % end
% % % 
% % % % Normalize common string containers
% % % if isstring(val)
% % %     v = cellstr(val);
% % % elseif ischar(val)
% % %     v = val;           % leave a single char row vector as-is
% % % elseif iscell(val)
% % %     % If it's a cell of strings/char, ensure cellstr form
% % %     if all(cellfun(@(x)ischar(x)||isstring(x), val))
% % %         v = cellfun(@char, val, 'UniformOutput', false);
% % %     else
% % %         v = val;       % heterogeneous cell: leave as-is
% % %     end
% % % else
% % %     v = val;           % numerics / logicals / etc.
% % % end
% % % end

