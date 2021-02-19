function smp = sampleHN(thetaD2,I,n)
if(~isnumeric(I))
    I = I.mu + I.Delta*abs(randn(n,1)) + sqrt(I.Gamma)*randn(n,1);
end
sigma=sigmaVec(thetaD2,I);
smp = I - sigma.*abs(randn(n,1));
end

