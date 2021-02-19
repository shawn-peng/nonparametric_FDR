s1=1:(10-1)/1000:10;
%s1=[repmat(1,1,500),repmat(10,1,500)];
s1=s1(:);
mu= s1;
sigmaSq = (s1).^2;
%sigmaSq = (s1);
%sigmaSq = 1;

x=mu + sqrt(sigmaSq).*randn(length(s1),1);

p= ones(length(s1),1);

L1= 0.5 * log (2*pi* sigmaSq) + 0.5 * ((x-mu).^2)./sigmaSq;
loss= mean(p .* L1);

%[netMu,netSigma]=HNNN(s1,x,p);
%NN=fitNN(s1,x,p);
[netMu,netSigma]=fitNNIter(s1,x,p);


XNew = reshape(s1', [1,1,size(x,2),size(x,1)]);

% Ymu=predict(netMu,XNew);
% Ysigma=predict(netSigma,XNew);
% muhat=Ymu(:,1);
% sigmahat=Ysigma(:,1);
muhat=netMu.predict(s1);
sigmahat=netSigma.predict(s1);
%[muhat,sigmahat]=NN.predict(s1);
subplot(1,2,1);
scatter(muhat,mu)
subplot(1,2,2)
scatter(sigmahat,sigmaSq)

