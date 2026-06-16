function o = overlay_opts(base, add)
o = base;
fn = fieldnames(add);
for i=1:numel(fn)
    o.(fn{i}) = add.(fn{i});
end
end
