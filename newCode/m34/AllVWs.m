function [Vc_1,Wc_1,Vc_2,Wc_2,VI1_1,WI1_1,VI1_2,WI1_2,VI2_2,WI2_2,VI2_3,WI2_3,VI3_3,WI3_3] = AllVWs(s1,s2,s3,zeta)

[Vc_1,Wc_1]=VW(s1,zeta.C);
[Vc_2,Wc_2]=VW(s2,zeta.C);
[VI1_1,WI1_1]=VW(s1,zeta.I1);
[VI1_2,WI1_2]=VW(s2,zeta.I1);
[VI2_2,WI2_2]=VW(s2,zeta.I2);
[VI2_3,WI2_3]=VW(s3,zeta.I2);
[VI3_3,WI3_3]=VW(s3,zeta.I3);
end

