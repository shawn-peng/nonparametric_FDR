function [p,p_c,p_I1] = s1Dens(x,zeta)
a=zeta.alpha;
p_c=a*snPdf(x,zeta.C);
p_I1=(1-a)*snPdf(x,zeta.I1);
p=p_c + p_I1;
end

