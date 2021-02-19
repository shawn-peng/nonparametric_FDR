function [p,p_c,p_I1,p_I2] = s2Dens(x,zeta)
a=zeta.alpha;
aC = cdfSN(zeta.C,x);
aI1I2=pBetweenI1I2(x,zeta);
p_c=a*aI1I2.*snPdf(x,zeta.C);
p_I1=a*(1-aC).*pdfI1(x,zeta);
p_I2=((1-a)+a*aC).*snPdf(x,zeta.I2);
p= p_c + p_I1 + p_I2;
end