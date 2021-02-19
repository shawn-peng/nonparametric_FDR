function [p1,p2,p] = pVec12(s1,zeta)
p1=snPdf(s1,zeta.C);
p2=snPdf(s1,zeta.I1);
a=zeta.alpha;
p=a*p1 + (1-a)*p2 ;
end

