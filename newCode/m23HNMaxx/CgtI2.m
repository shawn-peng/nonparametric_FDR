function p = CgtI2(thetaC,thetaI1,thetaD2)
n=10000;
xxC=SNSample(thetaC,n);
xxI1=SNSample(thetaI1,n);
sigma=  sigmaVec(thetaD2,xxI1);
xxI2 = xxI1-sigma.*randn(n,1);
p=sum(xxC>xxI2)/n;
end

