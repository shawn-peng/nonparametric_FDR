function [Vc_1,Wc_1,VI1_1,WI1_1] = AllVWs12HN(s1,zeta)

[Vc_1,Wc_1]=VW(s1,zeta.C);

[VI1_1,WI1_1]=VW(s1,zeta.I1);
end

