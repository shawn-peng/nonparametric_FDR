function [p,p_c,p_I1] = s1Dens_i(x,zeta)
a=zeta.alpha;

p_c=a*pdf(zeta.C,x);
p_I1=(1-a)*pdfI1(zeta,x);
%a2 = CgtI2(zeta.C,zeta.I1,zeta.D2);
p=p_c + p_I1;
end

