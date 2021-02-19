function s = topScore(zeta,n)

xC=SNSample(zeta.C,n);
xI1=SNSample(zeta.I1,n);
a=zeta.alpha;
n1=ceil(a*n);
s=[max(xI1(1:n1),xC(1:n1));xI1(n1+1:end)];
end

