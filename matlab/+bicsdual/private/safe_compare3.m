function [dR2,p] = safe_compare3(Xb,Xf,y)
% Nested OLS comparison on common finite rows; intercept included internally.
% Undefined fits/tests return NaN. No added rank returns DeltaR2=0, p=1.
dR2=NaN; p=NaN; y=y(:);
assert(size(Xb,1)==numel(y) && size(Xf,1)==numel(y), ...
    'bicsdual:Rows','Reduced/full/response row counts must agree.');
ok=isfinite(y)&all(isfinite(Xb),2)&all(isfinite(Xf),2);
y=y(ok); Xb=Xb(ok,:); Xf=Xf(ok,:); n=numel(y);
if n<2, return; end
Xb=[ones(n,1),Xb]; Xf=[ones(n,1),Xf];
TSS=sum((y-mean(y)).^2); if TSS==0, return; end
[~,eb,rb,Qb]=ols_projection(Xb,y);
[~,ef,rf,Qf]=ols_projection(Xf,y);
assert(norm(Qb-Qf*(Qf'*Qb),'fro')<=1e-7*max(1,sqrt(rb)), ...
    'bicsdual:NotNested','Reduced model is not nested in the full model at numerical tolerance.');
rss_b=sum(eb.^2); rss_f=sum(ef.^2);
gain=rss_b-rss_f;
assert(gain>=-1e-10*max(TSS,realmin), ...
    'bicsdual:Nonmonotonic','Full-model RSS exceeds reduced-model RSS.');
if rf==rb, dR2=0; p=1; return; end
dR2=max(0,min(1,gain/TSS));
df1=rf-rb; df2=n-rf;
if df1<=0 || df2<=0, return; end
if gain<=0, p=1; return; end
if rss_f==0, p=realmin; return; end
F=(gain/df1)/(rss_f/df2);
p=fcdf(F,df1,df2,'upper');
if p==0, p=realmin; end
end
