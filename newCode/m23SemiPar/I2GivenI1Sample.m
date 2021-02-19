function I2 = I2GivenI1Sample(zeta,I1)
n=10000;
I1=double(I1);
I2=nan(length(I1),1);
w=nan(length(I1),n);
x=random(zeta.I2,n);
xx=double(x);
[mu,ss]=predictMuAndSS(zeta.D1,x);
%p2=pdf(zeta.I2,x);

for i=1:n
    diff=I1-xx(i);
    jx=diff>=0;
    w(jx,i)= tnPdf(diff(jx),mu(i),sqrt(ss(i)));
    %w(jx,i)=w(jx,i)*p2(i);
    w(~jx,i)=0;
end

for j=1:length(I1)
   if sum(w(j,:))~=0
       I2(j) = xx(randsample(n,1,true,w(j,:)));
   else
       I2(j) = I1(j);
   end
end

I2=zeta.I2.bin(I2);

end

