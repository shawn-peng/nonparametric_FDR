function mat = simulateMaxHN(zeta,n)
%
alpha=zeta.alpha;
n1=int32(alpha*n);


xc = SNSample(zeta.C,n);
xI1 = SNSample(zeta.I1,n);
xI2 = sampleHN(zeta.D2, xI1, n); 
mx=max(xc(1:n1),xI1(1:n1));
mn=min(xc(1:n1),xI1(1:n1));
s1 = [mx;xI1(n1+1:n)];
s2 = [max(mn,xI2(1:n1));xI2(n1+1:n)];
mat=[s1,s2];
end
