function Q = Qfunc12(s1, zetaO, zetaN)
r=rVec12(s1,zetaO);
[Vc,Wc,VI1,WI1] = AllVWs12(s1,zetaO);

a=zetaO.alpha;

tmp1 = r*log(a) + (1-r)*log(1-a);
tmp2 = r.*qfun(s1,Vc,Wc,zetaN.C);
tmp3 =(1-r).*qfun(s1,VI1,WI1,zetaN.I1);
Q = sum(tmp1 + tmp2 + tmp3);

    function q= qfun(x,t,t2,theta)
        mu=theta.mu;Gamma=theta.Gamma;Delta=theta.Delta;
        q = log(Gamma) + ((x-mu).^2 - 2*Delta*(x-mu).*t + (Delta^2 + Gamma)*t2)/Gamma;
        q = -0.5*q;
    end
end

