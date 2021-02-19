function [zetaN,llN,i] = EM12(s1,sl1,sl2)
s1=s1(s1>0);
M=length(s1);
zeta0=initPar(s1);
zetaO=zeta0;
llO=-inf;

i=0;

while true
    QO=Qfunc12(s1,zetaO,zetaO);
    r=rVec12(s1,zetaO);
    [Vc,Wc,VI1,WI1] = AllVWs12(s1,zetaO);
    zetaN1 = zetaO;
    zetaN1.alpha = mean(r);
    zetaN=zetaN1;
    llN1 = loglikelihood12(s1,zetaN);
    QN1 = Qfunc12(s1,zetaO,zetaN);
    %if (llN1<llO||QN1<QO)
    if (QN1<QO)
        disp('non monotonic');
    end
   
    QN2=QN1;
    
    r_c = r;
    m_c = mfun(s1,Vc,zetaN.C.Delta); 
    zetaN3 = zetaN;
    zetaN3.C.mu=r_c'*m_c/sum(r_c);
    zetaN=zetaN3;
    llN3 = loglikelihood12(s1,zetaN);
    QN3 = Qfunc12(s1,zetaO,zetaN);
    %if (llN3<llN2||QN3<QN2)
    if (QN3<QN2)
        disp('non monotonic');
    end
   
    
    d_c = dfun(s1,Vc,zetaN.C.mu); 
    w_c = Wc;
    zetaN4 = zetaN;
    zetaN4.C.Delta=(r_c'*d_c)/(r_c'*w_c);
    zetaN=zetaN4;
    llN4 = loglikelihood12(s1,zetaN);
    QN4 = Qfunc12(s1,zetaO,zetaN);
    %if (llN4<llN3||QN4<QN3)
    if (QN4<QN3)
        disp('non monotonic');
    end
    
    g_c = gfun(s1,Vc,Wc,zetaN.C.mu,zetaN.C.Delta);
    zetaN5 = zetaN;
    zetaN5.C.Gamma=(r_c'*g_c)/sum(r_c);
    zetaN=zetaN5;
    llN5 = loglikelihood12(s1,zetaN);
    QN5 = Qfunc12(s1,zetaO,zetaN);
    %if (llN5<llN4||QN5<QN4)
    if (QN5<QN4)
        disp('non monotonic');
    end
    
    r_I1 = 1-r;
    m_I1 = mfun(s1,VI1,zetaN.I1.Delta); 
    zetaN6 = zetaN;
    zetaN6.I1.mu=r_I1'*m_I1/sum(r_I1);
    zetaN=zetaN6;
    llN6 = loglikelihood12(s1,zetaN);
    QN6 = Qfunc12(s1,zetaO,zetaN);
    %if (llN6<llN5||QN6<QN5)
    if (QN6<QN5)
        disp('non monotonic');
    end
   
    
    d_I1 = dfun(s1,VI1,zetaN.I1.mu); 
    w_I1 = WI1;
    zetaN7 = zetaN;
    zetaN7.I1.Delta=(r_I1'*d_I1)/(r_I1'*w_I1);
    zetaN=zetaN7;
    llN7 = loglikelihood12(s1,zetaN);
    QN7 = Qfunc12(s1,zetaO,zetaN);
    %if (llN7<llN6||QN7<QN6)
    if (QN7<QN6)
        disp('non monotonic');
    end
    
    g_I1 = gfun(s1,VI1,WI1,zetaN.I1.mu,zetaN.I1.Delta);
    zetaN8 = zetaN;
    zetaN8.I1.Gamma=(r_I1'*g_I1)/sum(r_I1);
    zetaN=zetaN8;
    llN = loglikelihood12(s1,zetaN);
    QN = Qfunc12(s1,zetaO,zetaN);
    %if (llN8<llN7||QN8<QN7)
    if (QN<QN7)
        disp('non monotonic');
    end
    
    if abs(llN-llO)<10^-7
        break 
    end
    llO=llN;
    %QO=QN;
    zetaO=zetaN;
    i=i+1;
end
    function m = mfun(x,v,Delta)
        m=x-v*Delta;
    end
    function d = dfun(x,v,mu)
        d=v.*(x-mu);
    end
    function g = gfun(x,v,w,mu,Delta)
        g=(x-mu).^2 - 2*v.*(x-mu).*Delta + w*Delta^2;
    end


    function zeta = initPar(s1)
        zeta=struct();
        zeta.alpha=0.5;
        zeta.C=momEst(s1(s1>median(s1)));
        zeta.C.Delta=sl1*abs(zeta.C.Delta);
        zeta.I1=momEst(s1(s1<=median(s1)));
         zeta.I1.Delta=sl2*abs(zeta.I1.Delta);
    end
end

