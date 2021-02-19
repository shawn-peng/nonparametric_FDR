function [r1,r2] = rVec34(s1,s2,s3,zeta)
[p1,p2,~,p]=pVec34(s1,s2,s3,zeta);
a=zeta.alpha;
b=zeta.beta;
r1=a*p1./p;
r1(isnan(r1))=0;
tmp=a*b*p1 + (1-a)*b*p2;
r2= (tmp)./p;
r2(tmp==0)=0;
end

