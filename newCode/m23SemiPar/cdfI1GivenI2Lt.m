function P = cdfI1GivenI2Lt(zeta,I1,up)
n=10000;
P=nan(length(I1),1);
I1=double(I1);
up=double(up);
xI2=random(zeta.I2,n);
xI2=double(xI2);
[mu,ss]=predictMuAndSS(zeta.D1,xI2);
xI1 = xI2 + TNSample(mu,sqrt(ss));

[xI1,ix]=sort(xI1);
xI2=xI2(ix);

for i=1:length(I1)
    ixx=searchLB(xI1,I1(i));
    P(i)=sum(xI2(1:ixx)<up(i))/n;
end

end
