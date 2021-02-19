function [p,p_c,p_I1,p_I2] = s2Dens(x,zeta)
a=zeta.alpha;
aI1 = cdfSN(zeta.I1,x);
aC = cdfSN(zeta.C,x);
aI2 = cdfSN(zeta.I2,x);
p_c=a*(1-aI1).*aI2.*snPdf(x,zeta.C);
p_I1=a*(1-aC).*snPdf(x,zeta.I1);
p_I2=((1-a)+a*aC).*snPdf(x,zeta.I2);
p= p_c + p_I1 + p_I2;
end