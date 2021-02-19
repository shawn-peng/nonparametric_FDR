function [zeta] = testHNEstimation(thetaD2,thetaI)
n=10000;
x1=SNSample(thetaI,n);
x2=I2Sample(thetaD2,x1,n);
[x1,ix]=sort(x1);
x2=x2(ix);
IX={};
J=1;
s_I2=x1-x2;
ss2=nan(length(x1),1);
for i=1:length(x1)
   ix=find(x1<x1(i)+J & x1>x1(i)-J) ;
   IX=[IX,ix];
   ss2(i)= sqrt(var([s_I2(ix);-s_I2(ix)])); 
end
net1=timedelaynet(0,5);
%net2=timedelaynet(0,5);
net1.trainParam.epochs=10;
Xi=cell(1,0);
Ai=cell(2,0);
zeta.D2.net=fitNN(x1,ss2,ones(length(ss2),1),net1);
x2Fit=I2Sample(zeta.D2,x1,n);
histogram(x2,'Normalization','pdf');
hold on
histogram(x2Fit, 'Normalization','pdf');

function net = fitNN(x,y,w,net)
       net= train(net,x',y',Xi,Ai,w');
    end
end

