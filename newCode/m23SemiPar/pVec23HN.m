function [p1,p2,p3,p,a1] = pVec23HN(s1,s2,zeta)
pc_1=zeta.C.pdf(s1);
pc_2=zeta.C.pdf(s2);
pI2_2=zeta.I2.pdf(s2);
pI1GI2_1=pdfI1GivenI2(zeta,s1,s2);
pI1AndI2Lt_1=pdfI1AndI2Lt(zeta,s1,s2);
pI1_2=pdfI1(zeta,s2);


a=zeta.alpha;
a1 = zeta.C.cdf(s2);
%Probability that I2 is greater than s2
%a2 = 1-tnCdf(s1-s2,mu,sqrt(ss));
p1=pc_1.*pI1_2;
p2=pI1AndI2Lt_1.*pc_2;
p3=pI1GI2_1.*pI2_2;
p= (1-a)*p3 +  a*(p1 + a1.*p3+ p2);
if any(isnan(p))||any(isinf(p))
    disp('hi')
end
end

