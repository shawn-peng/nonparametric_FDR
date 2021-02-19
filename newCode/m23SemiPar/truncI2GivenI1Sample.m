function I2 = truncI2GivenI1Sample(zeta,I1,u)
n=10000;
I1=double(I1);
u=double(u);
I2=nan(length(I1),1);
w=nan(length(I1),n);
x=random(zeta.I2,n);
xx=double(x);
[mu,ss]=predictMuAndSS(zeta.D1,x);
%p2=pdf(zeta.I2,x);

% I1_rep = repmat(I1, 1, n);
% u_rep = repmat(u, 1, n);
% xx_rep = repmat(xx', length(I1), 1);
% mu_rep = repmat(mu', length(I1), 1);
% ss_rep = repmat(ss', length(I1), 1);
% diff = I1_rep - xx_rep;
% flags = (u_rep - xx_rep) >= 0;
% 
% w = tnPdf(diff, mu_rep, sqrt(ss_rep));
% w(~flags) = 0;

for i=1:n
    diff=I1-xx(i);
    jx=diff>=I1-u;
    w(jx,i)= tnPdf(diff(jx),mu(i),sqrt(ss(i)));
    %w(jx,i)=w(jx,i)*p2(i);
    w(~jx,i)=0;
end

for j=1:length(I1)
   if sum(w(j,:))~=0
       I2(j) = xx(randsample(n,1,true,w(j,:)));
   else
       I2(j) = u(j);
   end
end

I2=zeta.I2.bin(I2);
end