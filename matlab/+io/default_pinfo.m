function pinfo = default_pinfo(N)
%DEFAULT_PINFO Minimal pinfo (parcellation) when no RSN/node metadata are supplied.
pinfo = struct();
pinfo.cis = {1:N};
pinfo.clabels_short = {'all'};
pinfo.nnode = N;
pinfo.rsn_id = ones(N,1);
pinfo.rsn = repmat({'all'},N,1);
end
