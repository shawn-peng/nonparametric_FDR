function [p,p_c,p_I1] = s1Dens(x,zeta)
a=zeta.alpha;
aI1 = cdfI1(x,zeta);
aC = cdfSN(zeta.C,x);
p_c=a*aI1.*snPdf(x,zeta.C);
p_I1=((1-a) + a*aC).*pdfI1(x,zeta);
%a2 = CgtI2(zeta.C,zeta.I1,zeta.D2);
p=p_c + p_I1;
end

