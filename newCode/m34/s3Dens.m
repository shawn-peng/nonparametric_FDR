function [p,p_I2,p_I3] = s3Dens(x,zeta)
a=zeta.alpha;
b=zeta.beta;
p_I2=(a+(1-a)*b)*snPdf(x,zeta.I2);
p_I3=(1-a)*(1-b)*snPdf(x,zeta.I3);
p= p_I2 + p_I3;
end