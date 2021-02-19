function [fdr,cdf,xm,t1p,fdr1p] = FDR(zeta)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
a=zeta.alpha;
n=10000;
xI=SNSample(zeta.I1,n);
xC=SNSample(zeta.C,n);
n1=ceil(a*n);
xm=[max(xI(1:n1),xC(1:n1));xI(n1+1:end)];
isI=xm==xI;
[xm,ix]=sort(xm);
isI=isI(ix);
cdf=(1:n)/n;
fdr=nan(length(xm),1);
for i=1:length(xm)
    fdr(i)=sum(isI(i:end))/n;
end
[~,im]=min(abs(fdr-0.01));
t1p=xm(im);
fdr1p=fdr(im);
end

