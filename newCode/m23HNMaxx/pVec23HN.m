function [p1,p2,p3,p,a1,a2] = pVec23HN(s1,s2,zeta)
pc_1=snPdf(s1,zeta.C);
pc_2=snPdf(s2,zeta.C);
pI1_1=snPdf(s1,zeta.I1);
pI1_2=snPdf(s2,zeta.I1);
[mu,ss]=predictMuAndSS(zeta.D2,s1);
pD2_12=tnPdf(s1-s2,mu,sqrt(ss));

a=zeta.alpha;
a1 = cdfSN(zeta.C,s2);
%Probability that I2 is greater than s2
a2 = 1-tnCdf(s1-s2,mu,sqrt(ss));
p1=pc_1.*pI1_2;
p2=pI1_1.*pc_2;
p3=pI1_1.*pD2_12;
p= (1-a)*p3 +  a*(p1 + a1.*p3+ a2.*p2);
if any(isnan(p))||any(isinf(p))
    disp('hi')
end
end

