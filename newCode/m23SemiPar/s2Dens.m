function [p,p_c,p_I1,p_I2] = s2Dens(x,zeta)
a=zeta.alpha;
aC = zeta.C.cdf(x);
aI1I2=pBetweenI1I2(x,zeta);
p_c=a*aI1I2.*zeta.C.pdf(x);
p_I1=a*(1-aC).*pdfI1(zeta,x);
p_I2=((1-a)+a*aC).*zeta.I2.pdf(x);
p= p_c + p_I1 + p_I2;
end