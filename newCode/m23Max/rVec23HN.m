function [r0,r11,r101,r100] = rVec23HN(s1,s2,zeta)
[p1,p2,p3,p]=pVec23HN(s1,s2,zeta);
a=zeta.alpha;
a1= cdfSN(zeta.C, s2);
a2=cdfSN(zeta.I2, s2);
r0=(1-a)*p3./p;
r11= a*p1./p;
r101=a*a2.*p2./p;
r100=a*a1.*p3./p;
r0(isnan(r0))=0;
r11(isnan(r11))= 0;
r101(isnan(r101))=0;
r100(isnan(r100))=0;
end

