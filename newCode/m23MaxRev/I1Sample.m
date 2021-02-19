function I1 = I1Sample(zeta,I2)
if isscalar(I2)
    I2=SNSample(zeta.I2,I2);
end
[mu,ss]=predictMuAndSS(zeta.D1,I2);
mu=zeros(length(I2),1); ss=ones(length(I2),1)+15;
delta=TNSample(mu,sqrt(ss));
I1=I2+delta;
end

