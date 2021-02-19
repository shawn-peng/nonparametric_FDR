function NN = fitNNSigma(x,y,p,mu,opts)
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
%nH=opts.H;
width=opts.W;
nH=length(width);
eta=opts.eta;

%SS=x^.2;

WIDTH=[1,width,1];
act1=@relu;
act2=@softPl;
aa={};aa(1:nH)={act1};
ACTS={aa{:}, act2};
y=y';
p=p';
mu=mu';


if isfield(opts,'NNSigma')
    par=opts.NNSigma.par;
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


    function [f,g]= loss(ss,ind)
        yk=y(ind);
        pk=p(ind);
        pk=pk./sum(pk);
        n=length(ss);
        muk = mu(ind);
        L1= 0.5 * log (2*pi* ss) + 0.5 * ((yk-muk).^2)./ss;
        f= sum(pk .* L1);
        gSS= 0.5./ss - 0.5 * ((yk-muk).^2)./(ss.^2);
        %gSS=(1/n)*gSS;
        g=pk.*gSS;
        if any(isnan(g),'all')
            disp('hi');
        end
    end
end

