function p = mg3Dens(x,thetaD3,thetaD2,thetaI)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
ss=10000;
p=nan(length(x),1);
z = thetaI.mu + thetaI.Delta*abs(randn(ss,1)) + sqrt(thetaI.Gamma)*randn(ss,1);
s2 = z - thetaD2.sigma*abs(randn(ss,1));
for i= 1: length(x)
    p(i)=mean(hnPdf(s2-x(i),thetaD3));
end
end