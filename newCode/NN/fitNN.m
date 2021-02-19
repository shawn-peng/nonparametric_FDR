function NN = fitNN(x,y,opts)
%DEF.H=3;
DEF.W=[2,3];
DEF.th=0.1;
DEF.EP=100;
DEF.eta=1;
if nargin > 2
opts=getOptions(opts,DEF);
else
    opts=DEF;
end

width=opts.W;
nH=length(width);
eta=opts.eta;

%SS=x^.2;

WIDTH=[1,width,1];
act1=@relu;
%act1=@sigmoid;
act2=@softPl;
aa={};aa(1:nH)={act1};
ACTS={aa{:}, act2};
y=y';

badInit=true;

lt=10;
[init,run,predict]=NNfuncsBN(WIDTH, ACTS, @loss);
while badInit
    i=1;
    if lt==100
        disp('error in NN initialization')
        break;
    end
while badInit
    if i>5
        break;
    end
par=init();
[NN.par,badInit]=run(par,x,opts.EP,eta,lt);
i=i+1;
end
lt=lt+5;
end

    
NN.predict=@(x)predict(NN.par,x);
NN.WIDTH=WIDTH;
NN.ACTS=ACTS;


    function [f,g]= loss(yhat,ind)
        yk=y(ind);
        n=length(yhat);
        L1= (yhat-yk).^2;
         f= mean(L1,2);
         g = (2/n)*(yhat-yk);
    end
end

