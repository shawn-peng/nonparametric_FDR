function NN = fitNNMu(x,y,p,ss,opts)
%DEF.H=3;
DEF.W=[2,3];
DEF.th=0.1;
DEF.EP=25;
DEF.eta=1;
if nargin > 4
opts=getOptions(opts,DEF);
else
    opts=DEF;
end

width=opts.W;
nH=length(width);
eta=opts.eta;

%SS=x^.2;

WIDTH=[1,width,1];
act1= @relu;
act2 = @softPl;
aa={};aa(1:nH)={act1};
ACTS={aa{:}, act2};
y=y';
p=p';
ss=ss';

if isfield(opts,'NNMu')
    par=opts.NNMu.par;
    ACTS=opts.NNMu.ACTS;
    WIDTH=opts.NNMu.WIDTH;
end

[init,run,predict]=NNfuncsBN(WIDTH, ACTS, @loss);

if ~exist('par','var')
    par=init();
end

NN.par=run(par,x,opts.EP,eta,10);
NN.predict=@(x)predict(NN.par,x);
NN.WIDTH=WIDTH;
NN.ACTS=ACTS;

    function [f,g]= loss(mu,ind)
        yk=y(ind);
        pk=p(ind);
        pk=pk./sum(pk);
        %n=length(mu);
        ssk = ss(ind);
        L1= 0.5 * log (2*pi* ssk) + 0.5 * ((yk-mu).^2)./ssk;
        f= sum(pk .* L1);
        gmu = -(yk-mu)./ssk;
        %gmu=(1/n)*gmu;
        g=pk.*gmu;
    end
end

