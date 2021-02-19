function P = pBetweenI1I2(x,zeta)
n=10000;
P=nan(length(x),1);
xI2=SNSample(zeta.I2,n);
[mu,ss]=predictMuAndSS(zeta.D1,xI2);
xI1 = xI2 + TNSample(mu,sqrt(ss));

[xI1,ix]=sort(xI1);
xI2=xI2(ix);

for i=1:length(x)
    ixx=searchUB(xI1,x(i));
    P(i)=sum(xI2(ixx:n)<x(i))/n;
end

end
