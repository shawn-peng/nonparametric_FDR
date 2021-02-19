function [Vc,Wc,VI1,WI1] = AllVWs(s1,zeta)
[Vc,Wc]=VW(s1,zeta.C);
[VI1,WI1]=VW(s1,zeta.I1);
end

