function [p,p_c,p_I1] = s1Dens(x,zeta)
a=zeta.alpha;
aI1 = cdfI1(zeta,x);
aC = cdf(zeta.C,x);
p_c=a*aI1.*pdf(zeta.C,x);
p_I1=((1-a) + a*aC).*pdfI1(zeta,x);
%a2 = CgtI2(zeta.C,zeta.I1,zeta.D2);
p=p_c + p_I1;
end

