function smp =I2SampleOld(thetaD2,I,n)
if(~isnumeric(I))
    I=SNSample(I,n);
end
sigma=sigmaVec(thetaD2,I);
smp=I-sigma.*abs(randn(n,1));
end