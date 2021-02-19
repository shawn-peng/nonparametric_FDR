function Q = Qfunc34(s1,s2,s3, zetaO, zetaN)
[r1,r2]=rVec34(s1,s2,s3,zetaO);
[Vc_1,Wc_1,Vc_2,Wc_2,VI1_1,WI1_1,VI1_2,WI1_2,VI2_2,WI2_2,VI2_3,WI2_3,VI3_3,WI3_3] = AllVWs(s1,s2,s3,zetaO);

a=zetaO.alpha;b=zetaO.beta;

tmp1 = r1*log(a) + (1-r1)*log(1-a) + r2*log(b) + (1-r2)*log(1-b);
tmp2 = r1.*qfun(s1,Vc_1,Wc_1,zetaN.C) + (1-r1).*r2.*qfun(s2,Vc_2,Wc_2,zetaN.C);
tmp3 = r1.*qfun(s2,VI1_2,WI1_2,zetaN.I1) + (1-r1).*qfun(s1,VI1_1,WI1_1,zetaN.I1);
tmp4 = (r1 + (1-r1).*r2).*qfun(s3,VI2_3,WI2_3,zetaN.I2) + (1-r1).*(1-r2).*qfun(s2,VI2_2,WI2_2,zetaN.I2);
tmp5 = (1-r1).*(1-r2).*qfun(s3,VI3_3,WI3_3,zetaN.I3);
Q = sum(tmp1 + tmp2 + tmp3 + tmp4 + tmp5);

    function q= qfun(x,t,t2,theta)
        mu=theta.mu;Gamma=theta.Gamma;Delta=theta.Delta;
        q = log(Gamma) + ((x-mu).^2 - 2*Delta*(x-mu).*t + (Delta^2 + Gamma)*t2)/Gamma;
        q = -0.5*q;
    end
end

