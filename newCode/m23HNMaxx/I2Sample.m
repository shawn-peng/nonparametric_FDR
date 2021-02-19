function I2 = I2Sample(zeta,I1)
if nargin < 2
    I1=SNSample(zeta.I1,n);
end
[mu,ss]=predictMuAndSS(zeta.D2,I1);
delta=TNSample(mu,sqrt(ss));
I2=I1-delta;
end

