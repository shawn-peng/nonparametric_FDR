function [p,p_c,p_I1,p_I2] = s2Dens(x,zeta)
a=zeta.alpha;
b=zeta.beta;
p_c=(1-a)*b*snPdf(x,zeta.C);
p_I1=a*snPdf(x,zeta.I1);
p_I2=(1-a)*(1-b)*snPdf(x,zeta.I2);
p= p_c + p_I1 + p_I2;
end