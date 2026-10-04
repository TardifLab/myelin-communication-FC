function test_bics_fdr_saved_results
%TEST_BICS_FDR_SAVED_RESULTS Synthetic integration checks; no real data needed.
% Run once after adding this folder to the MATLAB path.
% Checks known BH answers, ties, p=0/1/NaN, deduplication, family separation,
% RSN symmetry, masks, raw preservation, save/reload, and error guards.
root=tempname; mkdir(root); cleanup=onCleanup(@()rmdir(root,'s'));
w=warning('off','BICSFDR:ZeroP'); wc=onCleanup(@()warning(w));
p=[0;.001;.01;.025;.025;.5;.9;NaN];
expected=[0;.004;.08/3;.04;.04;2/3;1;NaN];
G=reshape(p,4,2); U=cell(1,8);
for m=1:8
    if m<=4, pg=G(m,:); else, pg=[1 1]; end
    U{m}=record(pg,m);
end
pairNames={'SPE-CMY','SPE-DE','NE-CMY','NE-DE'};
ids=[1 3 5;1 4 6;2 3 7;2 4 8]; cc=[1 5;1 6;2 5;2 6];
lev={'L1_route','L1_diff','L2_both'};
for x=1:4
    C=struct; C.stats_all.meta.CC=cc(x,:);
    for k=1:3
        r=U{ids(x,k)}; C.stats_all.(lev{k})=r.stats;
        C.P_all.(lev{k})=r.pblock; C.Delta_all.(lev{k})=r.dblock;
    end
    S=struct; S.meta.test_fixture=true; S.total=C; %#ok<STRNU>
    save(fullfile(root,['results_total_tmp' pairNames{x} '.mat']),'S');
end
F=bics_fdr_saved_results(root,'Save',false,'FCLabels',{'A','B'},'RSNLabels',{'X','Y'});
assert(height(F.family_summary)==6);
assert(isequal(F.family_summary.NTests,[8;24;32;8;24;32]));
assert(isequal(F.family_summary.NFDRSignificant,[5;15;20;0;0;0]));
Q=[F.pairs(1).conditions.total.stats_all.L1_route.global_q; ...
   F.pairs(3).conditions.total.stats_all.L1_route.global_q; ...
   F.pairs(1).conditions.total.stats_all.L1_diff.global_q; ...
   F.pairs(2).conditions.total.stats_all.L1_diff.global_q];
assert(isequal(isnan(Q(:)),isnan(expected)));
v=isfinite(expected); assert(max(abs(Q(v)-expected(v)))<1e-13);
assert(isequaln(F.pairs(1).conditions.total.stats_all.L1_route.global_q, ...
    F.pairs(2).conditions.total.stats_all.L1_route.global_q));
for x=1:4
    A=F.pairs(x).conditions.total;
    for k=1:3
        st=A.stats_all.(lev{k}); r=U{ids(x,k)};
        assert(isequaln(st.global_p,r.stats.global_p));
        assert(isequaln(st.node_deltaR2,r.stats.node_deltaR2));
        assert(isequaln(A.P_all.(lev{k}),r.pblock));
        for f=1:2
            q=squeeze(st.ntwk_q(1,:,:,f)); assert(isequaln(q,q'));
            assert(isequaln(reshape(A.Q_all.(lev{k})(:,1,f),[],1),[q(1,1);q(2,1);q(2,2)]));
            target=st.global_q(f);
            if isfinite(target)
                assert(max(abs(q(:)-target))<1e-13);
                nq=st.node_q(1,:,f); assert(max(abs(nq-target))<1e-13);
            else
                assert(all(isnan(q(:)))); assert(~any(st.node_h(1,:,f)));
            end
        end
    end
end
Fbefore=F;
[rawD,raw]=bics_fdr_view(F,2,'total','raw');
[maskedD,masked]=bics_fdr_view(F,2,'total','fdr');
assert(isequaln(F,Fbefore));
assert(all(isfinite(rawD.L2_both(:)))&&all(isnan(maskedD.L2_both(:))));
assert(all(isfinite(raw.L2_both.node_deltaR2(:)))&&all(isnan(masked.L2_both.node_deltaR2(:))));
assert(isequaln(raw.L2_both.node_p,masked.L2_both.node_p));
assert(strcmp(raw.meta.selected_condition,'total')&&strcmp(raw.meta.display_mode,'raw'));
assert(strcmp(masked.meta.selected_condition,'total')&&strcmp(masked.meta.display_mode,'fdr'));
% Repository exporter: unique models/RSN blocks, raw values, and q mapping.
paths=analysis.write_fdr_tables(F,fullfile(root,'tables'));
TG=readtable(paths{1,1}); TN=readtable(paths{1,2}); TV=readtable(paths{1,3});
assert(height(TG)==16 && height(TN)==48 && height(TV)==64);
assert(all(TN.rsn_j<=TN.rsn_i));
assert(sum(TG.reject_fdr)==5 && sum(TN.reject_fdr)==15 && sum(TV.reject_fdr)==20);
assert(all(isnan(TV.deltaR2_fdr(~TV.reject_fdr))));
keep = TV.reject_fdr == 1; % Logical mask, whether imported as numeric or logical
assert(isequal(TV.deltaR2_fdr(keep),TV.deltaR2(keep)));
for m=1:8
    select=strcmp(string(TG.communication_model),F.meta.model_order{m});
    assert(sum(select)==2);
end
must_error(@()analysis.write_fdr_tables(F,fullfile(root,'tables')),'BICSFDR:OutputExists');
Fsave=bics_fdr_saved_results(root,'OutputDir',fullfile(root,'out'));
z=load(Fsave.files.results,'F'); assert(isequaln(z.F,Fsave));
assert(exist(Fsave.files.summary,'file')==2&&exist(Fsave.files.global_tests,'file')==2);
must_error(@()bics_fdr_saved_results(root,'OutputDir',fullfile(root,'out')),'BICSFDR:OutputExists');
path=fullfile(root,'results_total_tmpSPE-DE.mat'); z=load(path,'S'); original=z.S;
S=original; S.total.stats_all.L1_route.global_p(1)=.002; save(path,'S');
must_error(@()bics_fdr_saved_results(root,'Save',false),'BICSFDR:DuplicateMismatch');
S=original; S.total.P_all.L1_route(1)=2; save(path,'S');
must_error(@()bics_fdr_saved_results(root,'Save',false),'BICSFDR:InvalidP');
S=original; S.total.stats_all.meta.CC=[2 6]; save(path,'S');
must_error(@()bics_fdr_saved_results(root,'Save',false),'BICSFDR:PairOrder');
S=original; S=rmfield(S,'total'); save(path,'S');
must_error(@()bics_fdr_saved_results(root,'Save',false),'BICSFDR:MissingCondition');
fprintf('All BICS FDR synthetic checks passed.\n');
clear wc cleanup
end

function r=record(pg,m)
Nfc=2; Nnode=4; Nrsn=2; map=[1 1;2 1;2 2];
st.global_p=pg; st.global_deltaR2=[.01 .02]+m*.001;
st.node_p=repmat(reshape(pg,[1 1 Nfc]),[1 Nnode 1]);
st.node_deltaR2=repmat(reshape(st.global_deltaR2,[1 1 Nfc]),[1 Nnode 1]);
st.ntwk_p=zeros(1,Nrsn,Nrsn,Nfc); st.ntwk_deltaR2=zeros(1,Nrsn,Nrsn,Nfc);
r.pblock=repmat(reshape(pg,[1 1 Nfc]),[3 1 1]);
r.dblock=repmat(reshape(st.global_deltaR2,[1 1 Nfc]),[3 1 1]);
for f=1:Nfc
    st.ntwk_p(1,:,:,f)=pg(f); st.ntwk_deltaR2(1,:,:,f)=st.global_deltaR2(f);
end
st.block_map=map; st.Nntwk=Nrsn; st.model_labels={'all_myelin'};
r.stats=st;
end

function must_error(fun,id)
try
    fun();
catch err
    assert(strcmp(err.identifier,id),'Expected %s, got %s: %s',id,err.identifier,err.message);
    return
end
error('Expected error %s but none occurred.',id);
end
