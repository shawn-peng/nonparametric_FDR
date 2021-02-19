function [p,P] = pdfAndCdfI2(thetaD2,thetaI,x)
n=10000;
z = SNSample(thetaI,n);
[mu,ss]=predictMuAndSS(thetaD2,z);

diff0=-mu./sqrt(ss);
P0=normcdf(diff0);

diff=-bsxfun(@minus,x,z');
ix=diff<0;
diff=diff-repmat(mu',length(x),1);
diff=diff./sqrt(ss');
pmat=normpdf(diff);
pmat(ix)=0;
pmat=pmat./(1-P0');
p=mean(pmat,2);

PMat=normcdf(diff);
PMat=PMat-repmat(P0',length(x),1);
PMat(ix)=0;
PMat=PMat./(1-P0');
P=1-mean(PMat,2);

end

