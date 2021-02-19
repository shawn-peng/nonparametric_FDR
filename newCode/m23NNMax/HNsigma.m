function sigma = HNSigma(s1,s2,eps)
[I1,ix1]=sort(s1);
%[I2,ix2]=sort(s2);
IX1={};
s_I2 = s1-s2; 
s_I2=s_I2(ix1);
J=0.1;
sigma=nan(length(I1),1);
for i = 1:length(I1)
   ix=find(I1<I1(i)+J & I1>I1(i)-J) ;
   IX1=[IX1,ix];
   sigma(i)= sqrt(var([s_I2(ix);-s_I2(ix)])); 
end
end

