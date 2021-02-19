function p = pdfI1GivenI2Lt(I1,up,zeta)
n=10000;
w=nan(length(I1),n);
x=SNSample(zeta.I2,n);
[mu,ss]=predictMuAndSS(zeta.D1,x);
%p2=snPdf(x,zeta.I2);


for i=1:n
    diff=I1-x(i);
    jx=diff>=I1-up;
    w(jx,i)= tnPdf(diff(jx),mu(i),sqrt(ss(i)));
    w(~jx,i)=0;
end
p=mean(w,2);
mn=min(p(p>0));
p(p==0)=mn;
end
