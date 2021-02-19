function p = I2Pdf(x,xx,theta)
%sigma=max(theta.a*xx + theta.b,10^-7);
[mu,ss]=predictMuAndSS(theta,xx);
p=tnPdf(x,mu,sqrt(ss));
end