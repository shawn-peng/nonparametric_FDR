function [zetaBest,llBest,lls,i] = EM12HNMax(mat,sl1,sl2,K,zeta0)
s1=mat(:,1);


M=length(s1);
if nargin < 5
zeta0=initPar(s1);
end

zetaO=zeta0;
llO=-inf;
%QO=-inf;
lls=[];
i=0;
llBest=-inf;

while true
    %QO=Qfunc23HN(s1,s2,zetaO,zetaO);
    [l0,l11,l10]=rVec12HN(s1,zetaO);
    [Vc_1,Wc_1,VI1_1,WI1_1] = AllVWs12HN(s1,zetaO);
    zetaN1 = zetaO;
    zetaN1.alpha = 1-mean(l0);
    zetaN=zetaN1;
    llN1 = loglikelihood12HN(s1,zetaN);
%     QN1 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN1<llO||QN1<QO)
%     if (QN1<QO)
%         disp('non monotonic');
%     end
  
    %xc1=truncAboveSNSample(zeta.C,s1);
    xc=SNSample(zetaO.C,M);
    xc2=truncAboveSNSample(zetaO.C,s1);
    [Vxc,Wxc]=VW(xc,zetaO.C);
    [Vxc_2,Wxc_2]=VW(xc2,zetaO.C);
    r_c = [l0; l11; l10];
    m_c = [mfun(xc,Vxc,zetaN.C.Delta); mfun(s1,Vc_1,zetaN.C.Delta); mfun(xc2,Vxc_2,zetaN.C.Delta)]; 
    zetaN3 = zetaN;
    zetaN3.C.mu=r_c'*m_c/sum(r_c);
    zetaN=zetaN3;
    llN3 = loglikelihood12HN(s1,zetaN);
%     QN3 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN3<llN2||QN3<QN2)
%     if (QN3<QN2)
%         disp('non monotonic');
%     end
   
    
    d_c = [dfun(xc,Vxc,zetaN.C.mu);dfun(s1,Vc_1,zetaN.C.mu);dfun(xc2,Vxc_2,zetaN.C.mu);]; 
    w_c = [Wxc;Wc_1;Wxc_2];
    zetaN4 = zetaN;
    zetaN4.C.Delta=(r_c'*d_c)/(r_c'*w_c);
    zetaN=zetaN4;
    llN4 = loglikelihood12HN(s1,zetaN);
%     QN4 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN4<llN3||QN4<QN3)
%     if (QN4<QN3)
%         disp('non monotonic');
%     end
    
    g_c = [gfun(xc,Vxc,Wxc,zetaN.C.mu,zetaN.C.Delta);gfun(s1,Vc_1,Wc_1,zetaN.C.mu,zetaN.C.Delta); ...
        gfun(xc2,Vxc_2,Wxc_2,zetaN.C.mu,zetaN.C.Delta)];
    zetaN5 = zetaN;
    zetaN5.C.Gamma=(r_c'*g_c)/sum(r_c);
    zetaN=zetaN5;
    llN5 = loglikelihood12HN(s1,zetaN);
    %QN5 = Qfunc23HN(s1,s2,zetaO,zetaN);
    %if (llN5<llN4||QN5<QN4)
%     if (QN5<QN4)
%         disp('non monotonic');
%     end
%     
    xI=truncAboveSNSample(zetaO.I1,s1);
    [VxI,WxI]=VW(xI,zetaO.I1);
    r_I1 = r_c;
    m_I1 = [mfun(s1,VI1_1,zetaN.I1.Delta); mfun(xI,VxI,zetaN.I1.Delta);mfun(s1,VI1_1,zetaN.I1.Delta)]; 
    zetaN6 = zetaN;
    zetaN6.I1.mu=r_I1'*m_I1/sum(r_I1);
    zetaN=zetaN6;
    llN6 = loglikelihood12HN(s1,zetaN);
%     QN6 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN6<llN5||QN6<QN5)
%     if (QN6<QN5)
%         disp('non monotonic');
%     end
   
    
    d_I1 = [dfun(s1,VI1_1,zetaN.I1.mu); dfun(xI,VxI,zetaN.I1.mu);dfun(s1,VI1_1,zetaN.I1.mu)]; 
    w_I1 = [WI1_1; WxI;WI1_1];
    zetaN7 = zetaN;
    zetaN7.I1.Delta=(r_I1'*d_I1)/(r_I1'*w_I1);
    zetaN=zetaN7;
    llN7 = loglikelihood12HN(s1,zetaN);
%     QN7 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN7<llN6||QN7<QN6)
%     if (QN7<QN6)
%         disp('non monotonic');
%     end
    
    g_I1 = [gfun(s1,VI1_1,WI1_1,zetaN.I1.mu,zetaN.I1.Delta); gfun(xI,VxI,WxI,zetaN.I1.mu,zetaN.I1.Delta); ...
        gfun(s1,VI1_1,WI1_1,zetaN.I1.mu,zetaN.I1.Delta)];
    zetaN8 = zetaN;
    zetaN8.I1.Gamma=(r_I1'*g_I1)/sum(r_I1);
    zetaN=zetaN8;
    llN = loglikelihood12HN(s1,zetaN);
%     QN8 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN8<llN7||QN8<QN7)
%     if (QN8<QN7)
%         disp('non monotonic');
%     end
    
    
    %plotFit(mat,zetaN)
    if(rem(i,1000)==0)
        disp('hi')
    end
    %disp(llN-llO);
    if abs(llN-llO) < 10^-7 || i==K
        break 
    end
    llO=llN;
    lls=[lls;llN];
    if llBest<llN
        llBest=llN;
        zetaBest=zetaN;
    end
    %QO=QN;
    zetaO=zetaN;
    i=i+1;
    %disp(i);
    
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
   
    
 
    function zeta = initPar(s1)
        zeta=struct();
        zeta.alpha=0.5;
        zeta.beta=0.5;
        zeta.C=momEst(s1(s1>median(s1)));
        zeta.C.Delta=sl1*abs(zeta.C.Delta);
        zeta.I1=momEst(s1(s1<=median(s1)));
        zeta.I1.Delta=sl2*abs(zeta.I1.Delta);     
    end
end

