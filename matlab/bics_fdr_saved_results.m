function F = bics_fdr_saved_results(svdir, varargin)
%BICS_FDR_SAVED_RESULTS BH-FDR for saved bicsdual Stage A results.
% F = bics_fdr_saved_results(svdir,'FCLabels',ylbl)
%
% Reads four MAT files containing S, in this order:
%   SPE-CMY, SPE-DE, NE-CMY, NE-DE.
% Options (name/value):
%   FilePrefix   : 'results_total_tmp' (use 'results_' for older files)
%   PairLabels   : {'SPE-CMY','SPE-DE','NE-CMY','NE-DE'}
%                  File suffixes; ORDER MUST REMAIN as above.
%   Conditions   : {} = all standard conditions present in the files
%   Alpha        : 0.05; significance uses adjusted p < Alpha
%   FCLabels     : {} = FC1, FC2, ...; provide labels in SAVED order
%   RSNLabels    : {} = RSN1, RSN2, ...; provide labels in SAVED order
%   OutputDir    : fullfile(svdir,'fdr_comparison')
%   Save         : true; saves F, family summary CSV, global test CSV
%   Overwrite    : false; existing outputs are protected
%
% Families are SEPARATE for individual and combined communication models,
% and SEPARATE for each condition and each spatial scale. All FC targets
% and all locations at that scale are included. For 8 targets / 7 RSNs /
% 400 nodes, each family contains 32 / 896 / 12800 tests respectively.
% Duplicate individual results across pair files must agree EXACTLY.
% This implementation supports M=1 (all selected myelin predictors jointly).
% It stops on multiple predictor-menu specifications rather than discarding
% or silently regrouping them. RSN input must contain unique undirected
% blocks (one triangle including the diagonal).
%
% F.pairs(xx).conditions.(condition) retains unfiltered Delta_all, P_all,
% and stats_all, adding Q_all (adjusted block p), H_all (block rejection),
% and global_q/global_h, node_q/node_h, ntwk_q/ntwk_h within stats_all.
% Raw inputs are NEVER overwritten. Stage B arrays are not copied into F;
% continue loading SI_all / elasticity_all / pert_all from the source S.
% See bics_fdr_view and bics_fdr_preview for raw / corrected display copies.
%
% NaN p-values occupy a family slot conservatively as p=1 internally;
% their output q stays NaN and rejection is false. Inf/out-of-range p stops.
% No Statistics/Bioinformatics Toolbox is required. No regression/spin tests
% are rerun. FDR assumes the supplied p-values are valid; it cannot repair
% observation dependence, numerical underflow, or an incorrect model fit.

ip = inputParser;
addRequired(ip,'svdir',@(x)ischar(x)||isstring(x));
addParameter(ip,'FilePrefix','results_total_tmp',@(x)ischar(x)||isstring(x));
addParameter(ip,'PairLabels',{'SPE-CMY','SPE-DE','NE-CMY','NE-DE'});
addParameter(ip,'Conditions',{});
addParameter(ip,'Alpha',0.05,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<1);
addParameter(ip,'FCLabels',{});
addParameter(ip,'RSNLabels',{});
addParameter(ip,'OutputDir','',@(x)ischar(x)||isstring(x));
addParameter(ip,'Save',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Overwrite',false,@(x)islogical(x)&&isscalar(x));
parse(ip,svdir,varargin{:}); o=ip.Results;
svdir=char(svdir); pairs=textcell(o.PairLabels);
assert(numel(pairs)==4 && numel(unique(pairs))==4, ...
    'BICSFDR:PairOrder','Provide four unique suffixes in SPE-CMY, SPE-DE, NE-CMY, NE-DE order.');
conditions=textcell(o.Conditions);
levels={'L1_route','L1_diff','L2_both'};
standard={'int_none','int_caliber','int_ed','int_both','total'};
sources=cell(1,4); provenance=struct([]);
for x=1:4
    path=fullfile(svdir,[char(o.FilePrefix) pairs{x} '.mat']);
    assert(exist(path,'file')==2,'BICSFDR:MissingFile','Missing file: %s',path);
    z=load(path,'S');
    assert(isfield(z,'S')&&isstruct(z.S)&&isscalar(z.S), ...
        'BICSFDR:MissingS','File must contain scalar struct S: %s',path);
    sources{x}=z.S;
    info=dir(path);
    provenance(x).path=path;
    provenance(x).bytes=info.bytes;
    provenance(x).modified=info.date;
    if isfield(z.S,'meta'), provenance(x).meta=z.S.meta; end
end
if isempty(conditions)
    for c=1:numel(standard)
        present=cellfun(@(s)isfield(s,standard{c}),sources);
        assert(all(present)||~any(present),'BICSFDR:MissingCondition', ...
            'Condition %s is present in only some files. Supply matching source files.',standard{c});
        if all(present), conditions{end+1}=standard{c}; end %#ok<AGROW>
    end
end
assert(~isempty(conditions)&&numel(unique(conditions))==numel(conditions), ...
    'BICSFDR:Conditions','No conditions found, or duplicate condition names supplied.');

F=struct;
F.meta.version='1.0.0';
F.meta.created=datestr(now,30);
F.meta.alpha=o.Alpha;
F.meta.method='Benjamini-Hochberg';
F.meta.family_rule='Separate condition x spatial scale x individual/combined; all FCs and locations.';
F.meta.nan_policy='NaN p occupies a family slot as p=1; output q=NaN and reject=false.';
F.meta.conditions=conditions;
F.meta.model_order={'SPE','NE','CMY','DE','SPE-CMY','SPE-DE','NE-CMY','NE-DE'};
F.meta.sources=provenance;
% Maps the three levels in each pair file to the eight unique configurations.
modelID=[1 3 5; 1 4 6; 2 3 7; 2 4 8];
expectedCC=[1 5;1 6;2 5;2 6];
refs=[1 1;3 1;1 2;2 2;1 3;2 3;3 3;4 3];
summaryRows={}; globalRows={}; zeroCount=0;
first=true;
for c=1:numel(conditions)
    cond=conditions{c}; R=cell(4,3);
    for x=1:4
        assert(isfield(sources{x},cond),'BICSFDR:MissingCondition', ...
            'Condition %s missing from pair %s.',cond,pairs{x});
        C=sources{x}.(cond);
        if isfield(C,'stats_all')&&isfield(C.stats_all,'meta')&&isfield(C.stats_all.meta,'CC')
            assert(isequal(C.stats_all.meta.CC(:)',expectedCC(x,:)), ...
                'BICSFDR:PairOrder','Saved CC does not match required pair order at %s / %s.',pairs{x},cond);
        end
        for k=1:3
            R{x,k}=read_record(C,levels{k},sprintf('%s / %s / %s',pairs{x},cond,levels{k}));
            r=R{x,k};
            if first
                Nfc=r.Nfc; Nnode=r.Nnode; Nrsn=r.Nrsn; B=r.B;
                blockMap=r.block_map; first=false;
                fc=textcell(o.FCLabels); rsn=textcell(o.RSNLabels);
                if isempty(fc), fc=arrayfun(@(i)sprintf('FC%d',i),1:Nfc,'UniformOutput',false); end
                if isempty(rsn), rsn=arrayfun(@(i)sprintf('RSN%d',i),1:Nrsn,'UniformOutput',false); end
                assert(numel(fc)==Nfc&&numel(rsn)==Nrsn,'BICSFDR:Labels', ...
                    'FCLabels/RSNLabels lengths do not match saved arrays.');
                F.meta.FCLabels=fc; F.meta.RSNLabels=rsn;
                F.meta.Nfc=Nfc; F.meta.Nnode=Nnode; F.meta.Nrsn=Nrsn;
                F.meta.block_map=blockMap;
            else
                assert(r.Nfc==Nfc&&r.Nnode==Nnode&&r.Nrsn==Nrsn&& ...
                    isequal(r.block_map,blockMap),'BICSFDR:LayoutMismatch', ...
                    'Dimensions/block ordering differ at %s / %s / %s.',pairs{x},cond,levels{k});
            end
        end
    end
    % Do not silently choose one copy if nominally identical tests differ.
    assert_duplicate(R{1,1},R{2,1},[cond ' SPE']);
    assert_duplicate(R{3,1},R{4,1},[cond ' NE']);
    assert_duplicate(R{1,2},R{3,2},[cond ' CMY']);
    assert_duplicate(R{2,2},R{4,2},[cond ' DE']);
    U=cell(1,8);
    for m=1:8, U{m}=R{refs(m,1),refs(m,2)}; end
    pg=zeros(8,Nfc); pb=zeros(8,B,Nfc); pn=zeros(8,Nnode,Nfc);
    for m=1:8
        pg(m,:)=U{m}.pg;
        pb(m,:,:)=reshape(U{m}.pb,[1 B Nfc]);
        pn(m,:,:)=U{m}.pn;
    end
    qg=nan(size(pg)); qb=nan(size(pb)); qn=nan(size(pn));
    groupNames={'individual','combined'};
    for g=1:2
        idx=(1:4)+(g-1)*4;
        pset={pg(idx,:),pb(idx,:,:),pn(idx,:,:)};
        qset=cell(1,3); scales={'global','RSN','node'};
        for s=1:3
            p=pset{s}; q=bh_adjust(p); qset{s}=q;
            h=isfinite(q)&q<o.Alpha;
            valid=isfinite(p); raw=valid&p<o.Alpha;
            zeroCount=zeroCount+sum(p(:)==0);
            pv=p(h); if isempty(pv), cutoff=NaN; else, cutoff=max(pv); end
            summaryRows(end+1,:)={cond,groupNames{g},scales{s},numel(p), ...
                sum(valid(:)),sum(~valid(:)),sum(raw(:)),sum(h(:)), ...
                sum(raw(:)&~h(:)),cutoff}; %#ok<AGROW>
        end
        qg(idx,:)=qset{1}; qb(idx,:,:)=qset{2}; qn(idx,:,:)=qset{3};
    end
    for m=1:8
        family=groupNames{1+(m>4)};
        for f=1:Nfc
            globalRows(end+1,:)={cond,family,F.meta.model_order{m},fc{f}, ...
                U{m}.dg(f),pg(m,f),qg(m,f),isfinite(qg(m,f))&&qg(m,f)<o.Alpha}; %#ok<AGROW>
        end
    end
    for x=1:4
        C=sources{x}.(cond);
        A=struct('Delta_all',C.Delta_all,'P_all',C.P_all,'stats_all',C.stats_all);
        A.Q_all=struct; A.H_all=struct;
        for k=1:3
            lev=levels{k}; m=modelID(x,k); st=A.stats_all.(lev);
            qblock=reshape(qb(m,:,:),[B 1 Nfc]);
            A.Q_all.(lev)=qblock;
            A.H_all.(lev)=isfinite(qblock)&qblock<o.Alpha;
            st.global_q=reshape(qg(m,:),size(st.global_p));
            st.global_h=isfinite(st.global_q)&st.global_q<o.Alpha;
            st.node_q=reshape(qn(m,:,:),size(st.node_p));
            st.node_h=isfinite(st.node_q)&st.node_q<o.Alpha;
            st.ntwk_q=nan(size(st.ntwk_p));
            for b=1:B
                i=blockMap(b,1); j=blockMap(b,2);
                st.ntwk_q(1,i,j,:)=reshape(qblock(b,1,:),[1 1 1 Nfc]);
                st.ntwk_q(1,j,i,:)=reshape(qblock(b,1,:),[1 1 1 Nfc]);
            end
            st.ntwk_h=isfinite(st.ntwk_q)&st.ntwk_q<o.Alpha;
            st.fdr_alpha=o.Alpha; st.fdr_method='BH';
            st.fdr_family=groupNames{1+(m>4)};
            st.fdr_n_global=4*Nfc; st.fdr_n_RSN=4*B*Nfc; st.fdr_n_node=4*Nnode*Nfc;
            A.stats_all.(lev)=st;
        end
        F.pairs(x).label=pairs{x};
        F.pairs(x).source_file=provenance(x).path;
        F.pairs(x).conditions.(cond)=A;
    end
    fprintf('BH-FDR complete: %s (each family: %d global, %d RSN, %d node tests).\n', ...
        cond,4*Nfc,4*B*Nfc,4*Nnode*Nfc);
end
F.family_summary=cell2table(summaryRows,'VariableNames', ...
    {'Condition','Family','Scale','NTests','NValid','NMissing','NRawSignificant', ...
     'NFDRSignificant','NLostAfterFDR','LargestRejectedRawP'});
F.global_tests=cell2table(globalRows,'VariableNames', ...
    {'Condition','Family','Model','FC','DeltaR2','RawP','AdjustedP','RejectFDR'});
F.meta.n_zero_p_unique=zeroCount;
if zeroCount>0
    warning('BICSFDR:ZeroP', ...
        ['%d stored p-values equal zero. They are retained; FDR cannot recover ' ...
         'precision lost by the original 1-fcdf calculation.'],zeroCount);
end
F.files=struct;
if o.Save
    outdir=char(o.OutputDir);
    if isempty(outdir), outdir=fullfile(svdir,'fdr_comparison'); end
    if ~exist(outdir,'dir'), mkdir(outdir); end
    F.files.results=fullfile(outdir,'BICS_FDR_results.mat');
    F.files.summary=fullfile(outdir,'BICS_FDR_family_summary.csv');
    F.files.global_tests=fullfile(outdir,'BICS_FDR_global_tests.csv');
    paths=struct2cell(F.files);
    for j=1:numel(paths)
        assert(o.Overwrite||exist(paths{j},'file')~=2,'BICSFDR:OutputExists', ...
            'Output exists: %s. Choose another OutputDir or set Overwrite=true.',paths{j});
        assert(~any(strcmp(paths{j},{provenance.path})),'BICSFDR:SourceOverwrite', ...
            'Refusing to overwrite a source file.');
    end
    tmp=[tempname(outdir) '.mat']; cleanup=onCleanup(@()delete_if_present(tmp));
    save(tmp,'F','-v7');
    [ok,msg]=movefile(tmp,F.files.results,'f'); assert(ok,'BICSFDR:Save','%s',msg);
    clear cleanup
    write_table_atomic(F.family_summary,F.files.summary,outdir);
    write_table_atomic(F.global_tests,F.files.global_tests,outdir);
    fprintf('Saved results to %s\n',F.files.results);
end
end

function r=read_record(C,lev,context)
assert(all(isfield(C,{'Delta_all','P_all','stats_all'})), ...
    'BICSFDR:MissingOutput','Missing Stage A outputs at %s.',context);
assert(isfield(C.Delta_all,lev)&&isfield(C.P_all,lev)&&isfield(C.stats_all,lev), ...
    'BICSFDR:MissingLevel','Missing %s.',context);
st=C.stats_all.(lev);
fields={'global_p','global_deltaR2','node_p','node_deltaR2', ...
    'ntwk_p','ntwk_deltaR2','block_map','Nntwk'};
assert(all(isfield(st,fields)),'BICSFDR:MissingStats','Incomplete stats at %s.',context);
assert(size(st.global_p,1)==1&&size(st.node_p,1)==1&&size(st.ntwk_p,1)==1&& ...
    size(C.P_all.(lev),2)==1,'BICSFDR:MultipleSpecs', ...
    'M must equal 1 at %s. Multiple predictor menus need an explicit family definition.',context);
r.Nfc=size(st.global_p,2); r.Nnode=size(st.node_p,2); r.Nrsn=st.Nntwk;
assert(isscalar(r.Nrsn)&&isfinite(r.Nrsn)&&r.Nrsn>=1&&r.Nrsn==fix(r.Nrsn), ...
    'BICSFDR:RSNCount','Invalid Nntwk at %s.',context);
r.block_map=st.block_map; r.B=size(r.block_map,1);
assert(isnumeric(r.block_map)&&size(r.block_map,2)==2&& ...
    all(isfinite(r.block_map(:)))&&all(r.block_map(:)>=1)&& ...
    all(r.block_map(:)<=r.Nrsn)&&all(r.block_map(:)==fix(r.block_map(:)))&& ...
    r.B==r.Nrsn*(r.Nrsn+1)/2&& ...
    size(unique(sort(r.block_map,2),'rows'),1)==r.B, ...
    'BICSFDR:BlockMap','Need unique undirected RSN blocks including diagonal at %s.',context);
r.pg=st.global_p; r.dg=st.global_deltaR2;
r.pn=st.node_p; r.dn=st.node_deltaR2;
r.pb=C.P_all.(lev); r.db=C.Delta_all.(lev);
require_shape(r.pg,[1 r.Nfc],context); require_shape(r.dg,[1 r.Nfc],context);
require_shape(r.pn,[1 r.Nnode r.Nfc],context); require_shape(r.dn,[1 r.Nnode r.Nfc],context);
require_shape(r.pb,[r.B 1 r.Nfc],context); require_shape(r.db,[r.B 1 r.Nfc],context);
require_shape(st.ntwk_p,[1 r.Nrsn r.Nrsn r.Nfc],context);
require_shape(st.ntwk_deltaR2,[1 r.Nrsn r.Nrsn r.Nfc],context);
ps={r.pg,r.pn,r.pb}; ds={r.dg,r.dn,r.db};
for k=1:3
    p=ps{k}; d=ds{k};
    assert(~any(isinf(p(:)))&&all(isnan(p(:))|(p(:)>=0&p(:)<=1)), ...
        'BICSFDR:InvalidP','p-values outside [0,1] or Inf at %s.',context);
    assert(~any(isinf(d(:)))&&all(isnan(d(:))|(d(:)>=-1e-12&d(:)<=1+1e-12)), ...
        'BICSFDR:InvalidDelta','Invalid DeltaR2 at %s.',context);
    assert(~any(isfinite(p(:))&~isfinite(d(:))),'BICSFDR:MissingDelta', ...
        'Finite p with missing DeltaR2 at %s.',context);
end
for b=1:r.B
    i=r.block_map(b,1); j=r.block_map(b,2);
    assert(isequaln(reshape(r.pb(b,1,:),[],1),reshape(st.ntwk_p(1,i,j,:),[],1))&& ...
        isequaln(reshape(r.pb(b,1,:),[],1),reshape(st.ntwk_p(1,j,i,:),[],1))&& ...
        isequaln(reshape(r.db(b,1,:),[],1),reshape(st.ntwk_deltaR2(1,i,j,:),[],1))&& ...
        isequaln(reshape(r.db(b,1,:),[],1),reshape(st.ntwk_deltaR2(1,j,i,:),[],1)), ...
        'BICSFDR:BlockCopies','Block vectors and symmetric stats matrices disagree at %s.',context);
end
if isfield(st,'model_labels'), r.model_labels=st.model_labels; else, r.model_labels={}; end
end

function require_shape(a,want,context)
actual=size(a); n=max(numel(actual),numel(want));
actual(end+1:n)=1; want(end+1:n)=1;
assert(isnumeric(a)&&isreal(a)&&isequal(actual,want), ...
    'BICSFDR:Shape','Unexpected array shape/type at %s.',context);
end

function assert_duplicate(a,b,context)
names={'pg','dg','pn','dn','pb','db','model_labels'};
for j=1:numel(names)
    assert(isequaln(a.(names{j}),b.(names{j})),'BICSFDR:DuplicateMismatch', ...
        ['Duplicate individual-model outputs disagree for %s (%s). ' ...
         'Do not pool until file versions/settings have been reconciled.'],context,names{j});
end
end

function q=bh_adjust(p)
shape=size(p); p=p(:); missing=isnan(p); p(missing)=1;
m=numel(p); [s,order]=sort(p,'ascend');
adjusted=s.*(m./(1:m)');
for j=m-1:-1:1, adjusted(j)=min(adjusted(j),adjusted(j+1)); end
adjusted=min(adjusted,1);
q=zeros(m,1); q(order)=adjusted; q(missing)=NaN; q=reshape(q,shape);
end

function c=textcell(x)
if isempty(x), c={}; return; end
if ischar(x), c={x}; elseif isstring(x), c=cellstr(x); else, c=x; end
assert(iscell(c),'BICSFDR:Text','Expected text or cell array of text.');
for j=1:numel(c)
    if isstring(c{j})&&isscalar(c{j}), c{j}=char(c{j}); end
    assert(ischar(c{j})&&isrow(c{j})&&~isempty(c{j}), ...
        'BICSFDR:Text','Labels and condition names must be nonempty text.');
end
c=c(:)';
end

function write_table_atomic(T,path,outdir)
tmp=[tempname(outdir) '.csv']; cleanup=onCleanup(@()delete_if_present(tmp));
writetable(T,tmp);
[ok,msg]=movefile(tmp,path,'f'); assert(ok,'BICSFDR:Save','%s',msg);
clear cleanup
end

function delete_if_present(path)
if exist(path,'file')==2, delete(path); end
end
