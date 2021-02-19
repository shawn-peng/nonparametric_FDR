function [zetaN,llN,i] = EM34HN(mat,sl1,sl2)
s1=mat(:,1);
s2=mat(:,2);
s3=mat(:,3);
M=length(s1);
zeta0=initPar(s1,s2,s3);
zetaO=zeta0;
llO=-inf;
QO=-inf;

i=0;

while true
    QO=Qfunc34HN(s1,s2,s3,zetaO,zetaO);
    [r1,r2]=rVec34HN(s1,s2,s3,zetaO);
    [Vc_1,Wc_1,Vc_2,Wc_2,VI1_1,WI1_1,VI1_2,WI1_2] = AllVWs34HN(s1,s2,zetaO);
    zetaN1 = zetaO;
    zetaN1.alpha = mean(r1);
    zetaN=zetaN1;
    llN1 = loglikelihood34HN(s1,s2,s3,zetaN);
    QN1 = Qfunc34HN(s1,s2,s3,zetaO,zetaN);
    %if (llN1<llO||QN1<QO)
    if (QN1<QO)
        disp('non monotonic');
    end
   
    
    zetaN2 = zetaN;
    zetaN2.beta=mean(r2);
     zetaN=zetaN2;
    llN2 = loglikelihood34HN(s1,s2,s3,zetaN);
    QN2 = Qfunc34HN(s1,s2,s3,zetaO,zetaN);
    %if (llN2<llN1||QN2<QN1)
    if (QN2<QN1)
        disp('non monotonic');
    end
   
    
    r_c = [r1; (1-r1).*r2];
    m_c = [mfun(s1,Vc_1,zetaN.C.Delta); mfun(s2,Vc_2,zetaN.C.Delta)]; 
    zetaN3 = zetaN;
    zetaN3.C.mu=r_c'*m_c/sum(r_c);
    zetaN=zetaN3;
    llN3 = loglikelihood34HN(s1,s2,s3,zetaN);
    QN3 = Qfunc34HN(s1,s2,s3,zetaO,zetaN);
    %if (llN3<llN2||QN3<QN2)
    if (QN3<QN2)
        disp('non monotonic');
    end
   
    
    d_c = [dfun(s1,Vc_1,zetaN.C.mu); dfun(s2,Vc_2,zetaN.C.mu)]; 
    w_c = [Wc_1; Wc_2];
    zetaN4 = zetaN;
    zetaN4.C.Delta=(r_c'*d_c)/(r_c'*w_c);
    zetaN=zetaN4;
    llN4 = loglikelihood34HN(s1,s2,s3,zetaN);
    QN4 = Qfunc34HN(s1,s2,s3,zetaO,zetaN);
    %if (llN4<llN3||QN4<QN3)
    if (QN4<QN3)
        disp('non monotonic');
    end
    
    g_c = [gfun(s1,Vc_1,Wc_1,zetaN.C.mu,zetaN.C.Delta); gfun(s2,Vc_2,Wc_2,zetaN.C.mu,zetaN.C.Delta)];
    zetaN5 = zetaN;
    zetaN5.C.Gamma=(r_c'*g_c)/sum(r_c);
    zetaN=zetaN5;
    llN5 = loglikelihood34HN(s1,s2,s3,zetaN);
    QN5 = Qfunc34HN(s1,s2,s3,zetaO,zetaN);
    %if (llN5<llN4||QN5<QN4)
    if (QN5<QN4)
        disp('non monotonic');
    end
    
    r_I1 = [1-r1; r1];
    m_I1 = [mfun(s1,VI1_1,zetaN.I1.Delta); mfun(s2,VI1_2,zetaN.I1.Delta)]; 
    zetaN6 = zetaN;
    zetaN6.I1.mu=r_I1'*m_I1/sum(r_I1);
    zetaN=zetaN6;
    llN6 = loglikelihood34HN(s1,s2,s3,zetaN);
    QN6 = Qfunc34HN(s1,s2,s3,zetaO,zetaN);
    %if (llN6<llN5||QN6<QN5)
    if (QN6<QN5)
        disp('non monotonic');
    end
   
    
    d_I1 = [dfun(s1,VI1_1,zetaN.I1.mu); dfun(s2,VI1_2,zetaN.I1.mu)]; 
    w_I1 = [WI1_1; WI1_2];
    zetaN7 = zetaN;
    zetaN7.I1.Delta=(r_I1'*d_I1)/(r_I1'*w_I1);
    zetaN=zetaN7;
    llN7 = loglikelihood34HN(s1,s2,s3,zetaN);
    QN7 = Qfunc34HN(s1,s2,s3,zetaO,zetaN);
    %if (llN7<llN6||QN7<QN6)
    if (QN7<QN6)
        disp('non monotonic');
    end
    
    g_I1 = [gfun(s1,VI1_1,WI1_1,zetaN.I1.mu,zetaN.I1.Delta); gfun(s2,VI1_2,WI1_2,zetaN.I1.mu,zetaN.I1.Delta)];
    zetaN8 = zetaN;
    zetaN8.I1.Gamma=(r_I1'*g_I1)/sum(r_I1);
    zetaN=zetaN8;
    llN8 = loglikelihood34HN(s1,s2,s3,zetaN);
    QN8 = Qfunc34HN(s1,s2,s3,zetaO,zetaN);
    %if (llN8<llN7||QN8<QN7)
    if (QN8<QN7)
        disp('non monotonic');
    end
    
    r_I2 = [(1-r1).*(1-r2); r1; (1-r1).*r2];
    s_I2 = [sfun(s1-s2); sfun(s2-s3);sfun(s1-s3)]; 
    zetaN9 = zetaN;
    zetaN9.D2.sigma=sqrt(r_I2'*s_I2/sum(r_I2));
    zetaN=zetaN9;
    llN9 = loglikelihood34HN(s1,s2,s3,zetaN);
    QN9 = Qfunc34HN(s1,s2,s3,zetaO,zetaN);
    %if (llN9<llN8||QN9<QN8)
    if (QN9<QN8)
        disp('non monotonic');
    end
   
    
    r_I3 = (1-r1).*(1-r2);
    s_I3 = sfun(s2-s3); 
    zetaN10 = zetaN;
    zetaN10.D3.sigma=sqrt(r_I3'*s_I3/sum(r_I3));
    zetaN=zetaN10;
    llN = loglikelihood34HN(s1,s2,s3,zetaN);
    QN = Qfunc34HN(s1,s2,s3,zetaO,zetaN);
    %if (llN12<llN11||QN12<QN11)
    if (QN<QN9)
        disp('non monotonic');
    end
    
    if(rem(i,1000)==0)
        disp('hi')
    end
    disp(llN-llO);
    if abs(llN-llO)<10^-7
        break 
    end
    llO=llN;
    %QO=QN;
    zetaO=zetaN;
    i=i+1;
    disp(i);
    
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

    function s = sfun(x)
        s=x.^2;
    end


    function zeta = initPar(s1,s2,s3)
        zeta=struct();
        zeta.alpha=0.5;
        zeta.beta=0.5;
        zeta.C=momEst(s1(s1>median(s1)));
        zeta.I1=momEst(s1(s1<=median(s1)));
        zeta.D2.sigma=sqrt(sum((s1-s2).^2))/M;
        zeta.D3.sigma=sqrt(sum((s2-s3).^2))/M;
    end
end

