function [Vc_1,Wc_1,Vc_2,Wc_2,VI1_1,WI1_1,VI1_2,WI1_2] = AllVWs34HN(s1,s2,zeta)

[Vc_1,Wc_1]=VW(s1,zeta.C);
[Vc_2,Wc_2]=VW(s2,zeta.C);
[VI1_1,WI1_1]=VW(s1,zeta.I1);
[VI1_2,WI1_2]=VW(s2,zeta.I1);

end

