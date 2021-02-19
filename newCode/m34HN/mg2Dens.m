function p = mg2Dens(x,thetaD,thetaI)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
ss=10000;
p=nan(length(x),1);
z = thetaI.mu + thetaI.Delta*abs(randn(ss,1)) + sqrt(thetaI.Gamma)*randn(ss,1);
for i= 1: length(x)
p(i)=mean(hnPdf(z-x(i),thetaD));
end
end

