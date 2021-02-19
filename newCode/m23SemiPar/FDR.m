function [fdr,cdf,x,t1p,fdr1p] = FDR(x,zeta)
a=zeta.alpha;
n=length(x);
x=sort(x);
xI=I1Sample(zeta,n);
xC=sample(zeta.C,n);
n1=ceil(a*n);
xm=[max(xI(1:n1),xC(1:n1));xI(n1+1:end)];
fdr=nan(n,1);
cdf=nan(n,1);
for i=1:n
    deno= (n+1-i)/n;
    %ix = xC>=x(i);
    %num=a*min(deno,sum(xC(ix)>=xI(ix))/n);
    ix = xI>=x(i);
    num = (1-a)*sum(ix)/n + a* sum(xI(ix)>=xC(ix))/n;
    fdr(i)=num/deno;
    cdf(i)=sum(xm<=x(i))/n;
end
[~,im]=min(abs(fdr-0.01));
t1p=x(im);
fdr1p=fdr(im);
x=zeta.I2.bin(x);
end

