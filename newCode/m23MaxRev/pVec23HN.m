function [p1,p2,p3,p,a1] = pVec23HN(s1,s2,zeta)
pc_1=snPdf(s1,zeta.C);
pc_2=snPdf(s2,zeta.C);
pI2_2=snPdf(s2,zeta.I2);
pI1GI2_1=pdfI1GivenI2(s1,s2,zeta);
pI1GI2Lt_1=pdfI1GivenI2Lt(s1,s2,zeta);
pI1_2=pdfI1(s2,zeta);


a=zeta.alpha;
a1 = cdfSN(zeta.C,s2);
%Probability that I2 is greater than s2
%a2 = 1-tnCdf(s1-s2,mu,sqrt(ss));
p1=pc_1.*pI1_2;
p2=pI1GI2Lt_1.*pc_2;
p3=pI1GI2_1.*pI2_2;
p= (1-a)*p3 +  a*(p1 + a1.*p3+ p2);
if any(isnan(p))||any(isinf(p))
    disp('hi')
end
end

