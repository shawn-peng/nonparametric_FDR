function [p,p_I2,p_I3] = s3Dens(x,zeta)
a=zeta.alpha;
b=zeta.beta;
p_I2=(a+(1-a)*b)*mg2Dens(x,zeta.D2,zeta.I1);
p_I3=(1-a)*(1-b)*mg3Dens(x,zeta.D3,zeta.D2,zeta.I1);
p= p_I2 + p_I3;
end
