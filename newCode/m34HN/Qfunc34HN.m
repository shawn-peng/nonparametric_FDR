function Q = Qfunc34HN(s1,s2,s3, zetaO, zetaN)
[r1,r2]=rVec34HN(s1,s2,s3,zetaO);
[Vc_1,Wc_1,Vc_2,Wc_2,VI1_1,WI1_1,VI1_2,WI1_2] = AllVWs34HN(s1,s2,zetaO);

a=zetaO.alpha;b=zetaO.beta;

tmp1 = r1*log(a) + (1-r1)*log(1-a) + r2*log(b) + (1-r2)*log(1-b);
tmp2 = r1.*qfun(s1,Vc_1,Wc_1,zetaN.C) + (1-r1).*r2.*qfun(s2,Vc_2,Wc_2,zetaN.C);
tmp3 = r1.*qfun(s2,VI1_2,WI1_2,zetaN.I1) + (1-r1).*qfun(s1,VI1_1,WI1_1,zetaN.I1);
tmp4 = r1.*qfun2(s2-s3,zetaN.D2)+ (1-r1).*r2.*qfun2(s1-s3,zetaN.D2) + (1-r1).*(1-r2).*qfun2(s1-s2,zetaN.D2);
tmp5 = (1-r1).*(1-r2).*qfun2(s2-s3,zetaN.D3);
Q = sum(tmp1 + tmp2 + tmp3 + tmp4 + tmp5);

    function q= qfun(x,t,t2,theta)
        mu=theta.mu;Gamma=theta.Gamma;Delta=theta.Delta;
        q = log(Gamma) + ((x-mu).^2 - 2*Delta*(x-mu).*t + (Delta^2 + Gamma)*t2)/Gamma;
        q = -0.5*q;
    end
    function q= qfun2(x,theta)
        q = 2*log(theta.sigma) + (x/theta.sigma).^2;
        q = -0.5*q;
    end
end

