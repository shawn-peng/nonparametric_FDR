function p = hnPdf(x,xx,theta)
%sigma=max(theta.a*xx + theta.b,10^-7);
sigma=sigmaVec(theta,xx);
p=2*normpdf(x,0,sigma);
p(x<0)=0;
end

