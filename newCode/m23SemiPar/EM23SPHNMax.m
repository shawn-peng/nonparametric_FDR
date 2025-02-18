function [zetaBest,llBest,lls,i] = EM23SPHNMax(mat,K,Z)
s1=mat(:,1);
s2=mat(:,2);
llBest=-inf;
%s3=mat(:,3);
net1=timedelaynet(0,5);
net2=timedelaynet(0,5);
net1.trainParam.epochs=10;
net2.trainParam.epochs=10;
Xi=cell(1,0);
Ai=cell(2,0);
EP=100;
net1 = network(1,2,[1;0],[1; 0],[0 0; 1 0],[0 1]);
net2 = network(1,2,[1;0],[1; 0],[0 0; 1 0],[0 1]);

%[I1,ix1]=sort(s1);

%IX1={};
diff = s1-s2; 
%s_I2=s_I2(ix1);
J=2;
diff_ss=nan(length(s2),1);
diff_mu=nan(length(s2),1);

% s2_mat = s2*ones(size(s2))';
% range_low = s2_mat' - J;
% range_high = s2_mat' + J;
% 
% kernel_diff = (diff*ones(size(s2))');
% kernel_diff = kernel_diff(s2_mat > range_low & s2_mat < range_high);
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
    [zeta0,S1,S2]=initPar(s1,s2);
else
    S1=zeta0.C.bin(s1);
    S2=zeta0.C.bin(s2);
end
zetaO=zeta0;


llO=-inf;
%QO=-inf;
lls=[];
i=0;

while true
    plotFit(mat,zeta0)
    pause(0.001)

    %QO=Qfunc23HN(s1,s2,zetaO,zetaO);
    [r0,r11,r101,r100]=rVec23HN(S1,S2,zetaO);
    
    %[Vc_1,Wc_1,Vc_2,Wc_2,VI2_1,WI1_1,VI1_2,WI1_2] = AllVWs23HN(s1,s2,zetaO);
    zetaN1 = zetaO;
    zetaN1.alpha = 1-mean(r0);
    zetaN=zetaN1;
%    llN1 = loglikelihood23HN(S1,S2,zetaN);
%     QN1 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN1<llO||QN1<QO)
%     if (QN1<QO)
%         disp('non monotonic');
%     end
  
    xc_r=zetaO.C.random(M);
    xc_rl=zetaO.C.ltRandom(S2);
    
    r_c = [r0; r11; r101; r100];
    C=concat(concat(concat(xc_r,S1),S2),xc_rl);
    zetaN3 = zetaN;
    zetaN3.C=zetaO.C.fixedCompsFit(C,r_c);
    zetaN=zetaN3;
%    llN3 = loglikelihood23HN(S1,S2,zetaN);
%     QN3 = Qfunc23HN(s1,s2,zetaO,zetaN);
%     %if (llN3<llN2||QN3<QN2)
%     if (QN3<QN2)
%         disp('non monotonic');
%     end
   
    
   
    xI2 = I2GivenI1Sample(zetaN,S2);
    xI2T = truncI2GivenI1Sample(zetaN,S1,S2);
    r_I2 = [r0+r100;r11;r101];
    I2=concat(concat(S2,xI2),xI2T);
    zetaN6 = zetaN;
    zetaN6.I2=zetaN.I2.fixedCompsFit(I2,r_I2);
    zetaN=zetaN6;
%    llN6 = loglikelihood23HN(S1,S2,zetaN);

    
    r_I1 = [r0+r100; r11; r101];
    [tmp_mu,tmp_ss]=predictMuAndSS(zetaN.D1,I2);
    diffC= [diff;S2-xI2;S1-xI2T];
    diffNTC=nonTruncNSample(tmp_mu,tmp_ss,diffC);
    zetaN9 = zetaN;
    zetaN9.D1.mu= sum(diffNTC.*r_I1)/sum(r_I1);
    zetaN9.D1.ss= sum((diffNTC.^2).*r_I1)/sum(r_I1) - zetaN9.D1.mu^2;
    
    zetaN=zetaN9;
    subplot(1,2,1)
    [mu,ss]=predictMuAndSS(zetaN.D1,S2);
        scatter(S2,diff_ss)
        hold on
        scatter(S2,ss);
        hold off
        subplot(1,2,2)
        scatter(S2,diff_mu)
        hold on
        scatter(S2,mu);
        hold off
  %  llN9 = loglikelihood23HN(S1,S2,zetaN)
    pause(0.1);

    zetaN10 = zetaN;
   
    zetaN=zetaN10;
    llN = loglikelihood23HN(S1,S2,zetaN)
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
    

 
    function [zeta,S1,S2] = initPar(s1,s2)
        zeta=struct();
        zeta.alpha=0.5;
        zeta.beta=0.5;
        [dC,dI2,S1,S2]=densEst_hist(s1,s2);
        zeta.C=dC;
        zeta.I2=dI2;
        
%         ixx=s1>median(s1);
%         opts.W=[2,3];
%         zeta.D1.netMu=fitNN(s2(~ixx),diff_mu(~ixx),opts);
%         zeta.D1.netSS=fitNN(s2(~ixx),diff_ss(~ixx),opts);
        
        zeta.D1.mu=0;
        zeta.D1.ss=1;
%         disp('parameters Initialized')
%         subplot(1,2,1)
%         scatter(s2,diff_ss)
%         hold on
%         scatter(s2,zeta.D1.netSS.predict(s2));
%         hold off
%         subplot(1,2,2)
%         scatter(s2,diff_mu)
%         hold on
%         scatter(s2,zeta.D1.netMu.predict(s2));
%         hold off
    end
end

