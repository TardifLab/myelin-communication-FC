function font(font_size, font_name)
%FONT Apply font settings to the current figure using built-in graphics.
%
% Compatibility helper for the original topology_spectra plotting code.

if nargin < 1 || isempty(font_size), font_size = 12; end
if nargin < 2 || isempty(font_name)
    font_name = get(groot,'DefaultAxesFontName');
end
fig = gcf;
ax = findall(fig,'Type','axes');
if ~isempty(ax), set(ax,'FontSize',font_size,'FontName',font_name); end
text_objects = findall(fig,'Type','text');
if ~isempty(text_objects)
    set(text_objects,'FontSize',font_size,'FontName',font_name);
end
end
