function p = CgtI1(thetaC,thetaI1)
n=10000;
xxC=SNSample(thetaC,n);
xxI1=SNSample(thetaI1,n);
p=sum(xxC>xxI1)/n;
end

