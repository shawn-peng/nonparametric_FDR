function  PlotCorr(zeta,s1,s2)
n=10000;
xC=SNSample(zeta.C,n);
xI1=SNSample(zeta.I1,n);
[mu,ss]=predictMuAndSS(zeta.D2,xI1);
a=zeta.alpha;
n1=ceil(a*n);
ss1=[max(xI1(1:n1),xC(1:n1));xI1(n1+1:end)];
s_I2= TNSample(mu, ss);
xI2=xI1-s_I2;

ss2=[max(xI2(1:n1),min(xI1(1:n1),xC(1:n1)));xI2(n1+1:end)];

subplot(2,3,1)
scatter(s2,s1-s2);
xlabel('s2')
ylabel('s1-s2');
title('real');
subplot(2,3,2)
scatter(s1,s1-s2);
xlabel('s1')
ylabel('s1-s2')
title('real');
subplot(2,3,3)
histogram(s1,'Normalization','pdf')
hold on
histogram(s2,'Normalization','pdf')
hold off
title('real');

subplot(2,3,4)
scatter(ss2,ss1-ss2);
xlabel('s2');
ylabel('s1-s2');
title('simulated');
subplot(2,3,5)
scatter(ss1,ss1-ss2);
xlabel('s1');
ylabel('s1-s2');
title('simulated');
subplot(2,3,6)
histogram(ss1,'Normalization','pdf')
hold on
histogram(ss2,'Normalization','pdf')
hold off
title('simulated');

figure;
subplot(2,2,1)
scatter(s2(s2<22),s1(s2<22)-s2(s2<22));
xlabel('s2')
ylabel('s1-s2');
title('real');
subplot(2,2,2)
scatter(s2,s1-s2);
xlabel('s2')
ylabel('s1-s2');
title('real');


subplot(2,2,3)
scatter(ss2,ss1-ss2);
xlabel('s2');
ylabel('s1-s2');
title('simulated');
subplot(2,2,4)
scatter(ss2,ss1-ss2);
xlabel('s2');
ylabel('s1-s2');
title('simulated');
xlim([0,60]);








end

