function [zetaBest,llBest,lls,i] = EM23HNMax(mat,slc,slI2,K,zeta0)
s1=mat(:,1);
s2=mat(:,2);
llBest=-inf;
%s3=mat(:,3);
net1=timedelaynet(0,5);
%net2=timedelaynet(0,5);
net1.trainParam.epochs=10;
%net2.trainParam.epochs=10;
Xi=cell(1,0);
Ai=cell(2,0);
EP=100;
%net1 = network(1,2,[1;0],[1; 0],[0 0; 1 0],[0 1]);
%net2 = network(1,2,[1;0],[1; 0],[0 0; 1 0],[0 1]);

%[I1,ix1]=sort(s1);

%IX1={};
diff = s1-s2; 
%s_I2=s_I2(ix1);
J=2;
diff_ss=nan(length(s2),1);
diff_mu=nan(length(s2),1);
for i = 1:length(s1)
   ix=find(s2<s2(i)+J & s2>s2(i)-J) ;
   %IX1=[IX1,ix];
   diff_ss(i)= var(diff(ix)); 
   diff_mu(i)= mean(diff(ix)); 
end

% IX2={};
% for i = 1:length(I2)
%    ix=find(I2<I2(i)+J & I2>I2(i)-J) ;
%    IX2=[IX2,ix];
% end


%s_I3 = s2-s3; 
%s_I3=s_I3(ix2);

M=length(s1);
if nargin<5
    zeta0=initPar(s1,s2);
end
zetaO=zeta0;
llO=-inf;
%QO=-inf;
lls=[];
i=0;

while true
    %QO=Qfunc23HN(s1,s2,zetaO,zetaO);
    [r0,r11,r101,r100]=rVec23HN(s1,s2,zetaO);
    
    %[Vc_1,Wc_1,Vc_2,Wc_2,VI2_1,WI1_1,VI1_2,WI1_2] = AllVWs23HN(s1,s2,zetaO);
    zetaN1 = zetaO;
    zetaN1.alpha = 1-mean(r0);
    zetaN=zetaN1;
    llN1 = loglikelihood23HN(s1,s2,zetaN);
%     QN1 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN1<llO||QN1<QO)
%     if (QN1<QO)
%         disp('non monotonic');
%     end
  
    %xc1=truncAboveSNSample(zeta.C,s1);
    [Vc_1,Wc_1]=VW(s1,zetaO.C);
    [Vc_2,Wc_2]=VW(s2,zetaO.C);
    xc=SNSample(zetaO.C,M);
    xc2=truncAboveSNSample(zetaO.C,s2);
    [Vxc,Wxc]=VW(xc,zetaO.C);
    [Vxc_2,Wxc_2]=VW(xc2,zetaO.C);
    r_c = [r0; r11; r101; r100];
    m_c = [mfun(xc,Vxc,zetaN.C.Delta); mfun(s1,Vc_1,zetaN.C.Delta); mfun(s2,Vc_2,zetaN.C.Delta); mfun(xc2,Vxc_2,zetaN.C.Delta)]; 
    zetaN3 = zetaN;
    zetaN3.C.mu=r_c'*m_c/sum(r_c);
    zetaN=zetaN3;
    llN3 = loglikelihood23HN(s1,s2,zetaN);
%     QN3 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN3<llN2||QN3<QN2)
%     if (QN3<QN2)
%         disp('non monotonic');
%     end
   
    
    d_c = [dfun(xc,Vxc,zetaN.C.mu);dfun(s1,Vc_1,zetaN.C.mu); dfun(s2,Vc_2,zetaN.C.mu);dfun(xc2,Vxc_2,zetaN.C.mu);]; 
    w_c = [Wxc;Wc_1; Wc_2;Wxc_2];
    zetaN4 = zetaN;
    zetaN4.C.Delta=(r_c'*d_c)/(r_c'*w_c);
    zetaN=zetaN4;
    llN4 = loglikelihood23HN(s1,s2,zetaN);
%     QN4 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN4<llN3||QN4<QN3)
%     if (QN4<QN3)
%         disp('non monotonic');
%     end
    
    g_c = [gfun(xc,Vxc,Wxc,zetaN.C.mu,zetaN.C.Delta);gfun(s1,Vc_1,Wc_1,zetaN.C.mu,zetaN.C.Delta); ...
        gfun(s2,Vc_2,Wc_2,zetaN.C.mu,zetaN.C.Delta);gfun(xc2,Vxc_2,Wxc_2,zetaN.C.mu,zetaN.C.Delta)];
    zetaN5 = zetaN;
    zetaN5.C.Gamma=(r_c'*g_c)/sum(r_c);
    zetaN=zetaN5;
    llN5 = loglikelihood23HN(s1,s2,zetaN);
    %QN5 = Qfunc23HN(s1,s2,zetaO,zetaN);
    %if (llN5<llN4||QN5<QN4)
%     if (QN5<QN4)
%         disp('non monotonic');
%     end
%   
    %xI2=SNSample(zetaN.I2,M);
    %xI2T=truncAboveSNSample(zetaN.I2,s2);
    [VI2_2,WI2_2]=VW(s2,zetaN.I2);
    xI2 = I2GivenI1Sample(zetaN,s2);
    xI2T = truncI2GivenI1Sample(zetaN,s1,s2);
    [VxI2,WxI2] = VW(xI2,zetaN.I2);
    [VxI2T,WxI2T] = VW(xI2T,zetaN.I2);
    r_I2 = [r0+r100;r11;r101];
    m_I2 = [mfun(s2,VI2_2,zetaN.I2.Delta);mfun(xI2,VxI2,zetaN.I2.Delta);mfun(xI2T,VxI2T,zetaN.I2.Delta)]; 
    zetaN6 = zetaN;
    zetaN6.I2.mu=r_I2'*m_I2/sum(r_I2);
    zetaN=zetaN6;
    llN6 = loglikelihood23HN(s1,s2,zetaN);
%     QN6 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN6<llN5||QN6<QN5)
%     if (QN6<QN5)
%         disp('non monotonic');
%     end
   
    
    d_I2 = [dfun(s2,VI2_2,zetaN.I2.mu); dfun(xI2,VxI2,zetaN.I2.mu);dfun(xI2T,VxI2T,zetaN.I2.mu)]; 
    w_I2 = [WI2_2; WxI2; WxI2T];
    zetaN7 = zetaN;
    zetaN7.I2.Delta=(r_I2'*d_I2)/(r_I2'*w_I2);
    zetaN=zetaN7;
    llN7 = loglikelihood23HN(s1,s2,zetaN);
%     QN7 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN7<llN6||QN7<QN6)
%     if (QN7<QN6)
%         disp('non monotonic');
%     end
    
    g_I2 = [gfun(s2,VI2_2,WI2_2,zetaN.I2.mu,zetaN.I2.Delta); gfun(xI2,VxI2,WxI2,zetaN.I2.mu,zetaN.I2.Delta);gfun(xI2T,VxI2T,WxI2T,zetaN.I2.mu,zetaN.I2.Delta) ];
    zetaN8 = zetaN;
    zetaN8.I2.Gamma=(r_I2'*g_I2)/sum(r_I2);
    zetaN=zetaN8;
    llN8 = loglikelihood23HN(s1,s2,zetaN);
%     QN8 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN8<llN7||QN8<QN7)
%     if (QN8<QN7)
%         disp('non monotonic');
%     end
    
    r_I1 = [r0+r100; r11; r101];
    %r_I2 = [r0+r100;];
    opts.NNSigma=zetaN.D1.netSS;
    opts.NNMu=zetaN.D1.netMu;
    I2C=[s2;xI2;xI2T];
    tmp_ss=opts.NNSigma.predict(I2C);
    tmp_mu=opts.NNMu.predict(I2C);
    diffC= [diff;s2-xI2;s1-xI2T];
    diffNTC=nonTruncNSample(tmp_mu,tmp_ss,diffC);
    %s_I2C = [s_I2NT];
    zetaN9 = zetaN;
    zetaN9.D1.netMu=fitNNMu(I2C,diffNTC,r_I1, tmp_ss, opts);
    zetaN9.D1.netSS=fitNNSigma(I2C, diffNTC, r_I1, tmp_mu, opts);
    
    zetaN=zetaN9;
    subplot(1,2,1)
        scatter(s2,diff_ss)
        hold on
        scatter(s2,zetaN.D1.netSS.predict(s2));
        hold off
        subplot(1,2,2)
        scatter(s2,diff_mu)
        hold on
        scatter(s2,zetaN.D1.netMu.predict(2));
        hold off
    llN9 = loglikelihood23HN(s1,s2,zetaN)

    zetaN10 = zetaN;
   
    zetaN=zetaN10;
    llN = loglikelihood23HN(s1,s2,zetaN);
%     QN = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN12<llN11||QN12<QN11)
%     if (QN<QN9)
%         disp('non monotonic');
%     end
    %plotFit(mat,zetaN)
    if(rem(i,1000)==0)
        disp('hi')
    end
   % disp(llN-llO);
    if abs(llN-llO)<10^-7 || i==K
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
   % disp(i);
    
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
        ixx=s1>median(s1);
        zeta.C=momEst(s1(ixx));
        zeta.C.Delta = slc*abs(zeta.C.Delta);
        zeta.I2=momEst(s2(~ixx));
        zeta.I2.Delta = slI2*abs(zeta.I2.Delta);
%         ss2=nan(length(I1),1);
%         for jj = 1 : length(I1)
%             ss2(jj)= sqrt(var(s_I2(IX1{jj}))); 
%         end
        %[a,b]=linearFit(I1,ss2,ones(length(ss2),1));
        %zeta.D2.a=a;zeta.D2.b=b;
        opts.W=[2,3];
        zeta.D1.netMu=fitNN(s2(~ixx),diff_mu(~ixx),opts);
        zeta.D1.netSS=fitNN(s2(~ixx),diff_ss(~ixx),opts);
        disp('parameters Initialized')
        subplot(1,2,1)
        scatter(s2,diff_ss)
        hold on
        scatter(s2,zeta.D1.netSS.predict(s2));
        hold off
        subplot(1,2,2)
        scatter(s2,diff_mu)
        hold on
        scatter(s2,zeta.D1.netMu.predict(s2));
        hold off
%         ss3=nan(length(I2),1);
%         for jj = 1 : length(I2)
%             ss3(jj)= sqrt(var(s_I3(IX2{jj}))); 
%         end
%         %[a,b]=linearFit(I2,ss3,ones(length(ss3),1));
%         %zeta.D3.a=a;zeta.D3.b=b;
%         zeta.D3.net=fitNN(I2,ss3,ones(length(ss3),1),net2);
    end
end

