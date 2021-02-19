function r = rVec12(s1,zeta)
[p1,~,p]=pVec12(s1,zeta);
a=zeta.alpha;
r=a*p1./p;
r(isnan(r))=0;
end

