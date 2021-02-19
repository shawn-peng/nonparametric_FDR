function [r0,r11,r10] = rVec12HN(s1,zeta)
[p1,p2,p]=pVec12HN(s1,zeta);
a=zeta.alpha;
a1 = cdfSN(zeta.C, s1);
a2 = cdfSN(zeta.I1, s1);
r0 = (1-a)*p2./p;
r11 = a*a2.*p1./p;
r10 = a*a1.*p2./p;
r0(isnan(r0))=0;
r11(isnan(r11))= 0;
r10(isnan(r10))=0;
end

