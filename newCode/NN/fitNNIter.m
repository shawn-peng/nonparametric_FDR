function [NNMu,NNSigma] = fitNNIter(x,y,p,opts)
niter=2;
ss=ones(length(x),1);
mu=x;
ss=x.^2;
%ss=x;
scatter(x,y)
%mu=nan(length(x),1);
%opts.H=3;
opts.W=[2,3,5];
opts.eta=0.01;
%opts=struct();
opts.NNMu=fitNN(x,mu,opts);
opts.NNSigma=fitNN(x,ss,opts);
for i=1:niter
    if i > 1
        opts.NNMu=NNMu;
        opts.NNSigma=NNSigma;
    end
    NNMu=fitNNMu(x,y,p,ss,opts);
    mu=NNMu.predict(x);
    NNSigma=fitNNSigma(x,y,p,mu,opts);
    ss=NNSigma.predict(x);
end

end

