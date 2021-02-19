function [p1,p2,p3,p] = pVec34(s1,s2,s3,zeta)
% pc_1=snPdf(s1,zeta.C.mu,zeta.C.gamma,zeta.C.Del,true);
% pc_2=snPdf(s2,zeta.C.mu,zeta.C.gamma,zeta.C.Del,true);
% pI1_1=snPdf(s1,zeta.I1.mu,zeta.I1.gamma,zeta.I1.Del,true);
% pI1_2=snPdf(s2,zeta.I1.mu,zeta.I1.gamma,zeta.I1.Del,true);
% pI2_2=snPdf(s2,zeta.I2.mu,zeta.I2.gamma,zeta.I2.Del,true);
% pI2_3=snPdf(s3,zeta.I2.mu,zeta.I2.gamma,zeta.I2.Del,true);
% pI3_3=snPdf(s3,zeta.I3.mu,zeta.I3.gamma,zeta.I3.Del,true);
pc_1=snPdf(s1,zeta.C);
pc_2=snPdf(s2,zeta.C);
pI1_1=snPdf(s1,zeta.I1);
pI1_2=snPdf(s2,zeta.I1);
pI2_2=snPdf(s2,zeta.I2);
pI2_3=snPdf(s3,zeta.I2);
pI3_3=snPdf(s3,zeta.I3);
a=zeta.alpha;
b=zeta.beta;
p1=pc_1.*pI1_2.*pI2_3;
p2=pI1_1.*pc_2.*pI2_3;
p3=pI1_1.*pI2_2.*pI3_3;
p=a*p1 + (1-a)*b*p2 + (1-a)*(1-b)*p3;
end

