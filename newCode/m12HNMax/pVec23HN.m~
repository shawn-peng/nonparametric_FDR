function [p1,p2,p] = pVec12HN(s1,zeta)
% pc_1=snPdf(s1,zeta.C.mu,zeta.C.gamma,zeta.C.Del,true);
% pc_2=snPdf(s2,zeta.C.mu,zeta.C.gamma,zeta.C.Del,true);
% pI1_1=snPdf(s1,zeta.I1.mu,zeta.I1.gamma,zeta.I1.Del,true);
% pI1_2=snPdf(s2,zeta.I1.mu,zeta.I1.gamma,zeta.I1.Del,true);
% pI2_2=snPdf(s2,zeta.I2.mu,zeta.I2.gamma,zeta.I2.Del,true);
% pI2_3=snPdf(s3,zeta.I2.mu,zeta.I2.gamma,zeta.I2.Del,true);
% pI3_3=snPdf(s3,zeta.I3.mu,zeta.I3.gamma,zeta.I3.Del,true);
%n=10000;
%xC=SNSample(zeta.C,n);
%xI1=SNSample(zeta.I1,n);
%sigma=sigmaVec(zeta.D2,xI1);
%xI2=xI1- sigma.* randn(n,1);
pc_1=snPdf(s1,zeta.C);
pI1_1=snPdf(s1,zeta.I1);

a=zeta.alpha;
a1 = cdfSN(zeta.C,s1);
a2 = cdfSN(zeta.I1,s1);

p1=pc_1;
p2=pI1_1;

p= (1-a)*p2 +  a*(a2.*p1 + a1.*p2);
if any(isnan(p))||any(isinf(p))
    disp('hi')
end
end

