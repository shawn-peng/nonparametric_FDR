function[mat,s1,s2] = dataGen(zeta,n)
xC=SNSample(zeta.C,n);
xI1=SNSample(zeta.I1,n);
xI2=I2Sample(zeta,xI1);
a=zeta.alpha;
n1=ceil(a*n);
s1=[max(xI1(1:n1),xC(1:n1));xI1(n1+1:end)];
s2=[max(xI2(1:n1),min(xI1(1:n1),xC(1:n1)));xI2(n1+1:end)];
mat(:,1)=s1;
mat(:,2)=s2;
end

