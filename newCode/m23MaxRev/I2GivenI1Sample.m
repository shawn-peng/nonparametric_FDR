function I2 = I2GivenI1Sample(zeta,I1)
n=10000;
I2=nan(length(I1),1);
w=nan(length(I1),n);
x=SNSample(zeta.I2,n);
[mu,ss]=predictMuAndSS(zeta.D1,x);
p2=snPdf(x,zeta.I2);

for i=1:n
    diff=I1-x(i);
    jx=diff>=0;
    w(jx,i)= tnPdf(diff(jx),mu(i),sqrt(ss(i)));
    w(jx,i)=w(jx,i)*p2(i);
    w(~jx,i)=0;
end

for j=1:length(I1)
   if sum(w(j,:))~=0
       I2(j) = x(randsample(n,1,true,w(j,:)));
   else
       I2(j) = I1(j) - TNSample(0,2);
   end
end

end

