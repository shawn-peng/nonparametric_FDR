function NN = fitNNSigma4(x,y,p,MU,opts)
DEF.H=4;
DEF.W=10;
DEF.th=0.1;
DEF.EP=1000;
DEF.eta=1;
if nargin > 4
opts=getOptions(opts,DEF);
else
    opts=DEF;
end
nH=opts.H;
width=opts.W;
eta=opts.eta;

%SS=x^.2;

WIDTH=[1,width,width,width,width,1];
ACT={@act1, @act1, @act1, @act1,@act2};
X=[x';ones(1,length(x))];
y=y';
p=p';
MU=MU';

if isfield(opts,'NN_Sigma')
    Ws=opts.NN_Sigma.Ws;
else
    for i=1:(length(WIDTH)-1)
        Ws{i}=sqrt(2/WIDTH(i))*randn(WIDTH(i+1),WIDTH(i)+1);
    end
end

%Ws{1}=[2,-10;-1/4,5/4];
%Ws{2}=[1,0,0;0,1,0];
%Ws{3}=[1,1,0];
%Ws{1}=[-1/9,10/9;2,-10;];
%Ws{2}=[1,0,0;0,1,0];
%Ws{3}=[1,1,0];
losses=[];
for k=1:opts.EP
ind=getIX(length(x));
Xk=X(:,ind);
yk=y(ind);
pk=p(ind);
muk=MU(ind);

%IN={Xk};
F={};
G={};
FW={};
GW={};
IN{1}=Xk;

for i=1:(length(WIDTH)-1)
    %[f,g]=linear(W,f);
    %FW{i}=f;
    %GW{i}=g;
    f=Ws{i}*IN{i};
    act=ACT{i};
    [f,g]=act(f);
    F{i}=f;
    if any(isnan(f),'all')
        disp('hi');
    end
%     if(i<=2)
%         subplot(3,2,2*(i-1)+1)
%         scatter(x,f(1,:))
%         hold on;
%         subplot(3,2,2*i)
%         scatter(x,f(2,:))
%         hold on;
%     else
%         subplot(3,2,5)
%         scatter(x,f)
%         hold on;
%     end
%     if k==1
%         subplot(3,2,6)
%         histogram(f,'Normalization','pdf');
%     end
    G{i}=g;
    if i<(length(WIDTH)-1)
        IN{i+1}=[f;ones(1,size(f,2))];
        %G{i}=[g;zeros(1,size(g,2))];
    else
        IN{i+1}=f;
    end
end
WW=Ws;
% subplot(2,2,3)
% scatter(x,f)
% hold on;
[f,g]=gaussianLoss(f,yk,pk,muk);
disp(f);
losses=[losses,f];
if k>100
    if losses(k-100)-losses(k)<10^-3
        eta = max(eta/2,0.001);
    end
end
 if any(isnan(g),'all')
        disp('hi');
 end
del5=g.*G{5};
WW{5}=outProdSum(del5,IN{5});

% WW=Ws{i};
% %G=GW{i};
% for j = 1:size(WW,1)
% WW(j,:)= sum(t1(j,:).*squeeze(G(j,:,:)),2)'; 
% end
% Ws{3}= Ws{3} + eta*WW;
w=Ws{5};
w=w(:,1:(size(Ws{5},2)-1));
del4 = (w'*del5).* G{4};
WW{4}=outProdSum(del4,IN{4});
%Ws{2}= Ws{2} + eta*WW;

w=Ws{4};
w=w(:,1:(size(Ws{4},2)-1));
del3 = (w'*del4).* G{3};
WW{3}=outProdSum(del3,IN{3});
%Ws{2}= Ws{2} + eta*WW;

w=Ws{3};
w=w(:,1:(size(Ws{3},2)-1));
del2 = (w'*del3).* G{2};
WW{2}=outProdSum(del2,IN{2});

w=Ws{2};
w=w(:,1:(size(Ws{2},2)-1));
del1 = (w'*del2).* G{1};
WW{1}=outProdSum(del1,IN{1});

Ws{1}= Ws{1} - eta*WW{1};
Ws{2}= Ws{2} - eta*WW{2};
Ws{3}= Ws{3} - eta*WW{3};
Ws{4}= Ws{4} - eta*WW{4};
Ws{5}= Ws{5} - eta*WW{5};

end
NN.Ws=Ws;
NN.predict=@predict;

 function [ss]=predict(x)
        XX=[x';ones(1,length(x))];
        INP{1}=XX;
        for ii=1:(length(WIDTH)-1)
            %[f,g]=linear(W,f);
            %FW{i}=f;
            %GW{i}=g;
            f=Ws{ii}*INP{ii};
            act=ACT{ii};
            [f,~]=act(f);
            if ii<(length(WIDTH)-1)
                INP{ii+1}=[f;ones(1,size(f,2))];
                %G{i}=[g;zeros(1,size(g,2))];
            else
                INP{ii+1}=f;
            end
        end
        ss=f(1,:)';
    end

 
    function ops = outProdSum(v1,v2)
         ops = zeros(size(v1,1),size(v2,1));
        n=size(v1,2);
        for ix=1:n
            ops = ops + v1(:,ix)*v2(:,ix)';
        end
        ops=ops/n;
    end

%     function [f,g]=iden(X)
%         f=X;
%         g=ones(size(X));
%     end
%     function [f,g]=act1(X)
%         f=1./(1+exp(-X));
%         g=f.*(1-f);
%     end
 function [f,g]=act1(X)
     f=X;
     ix=X<0;  
     f(ix)= X(ix)*0.01;
     f(ix)=0;
     g=ones(size(X));
     g(ix)=0;
    end
%     function [f,g]=act2(X)
%         C=opts.th;
%         f=max(C,X);
%         ix=X<C;
%         g=ones(size(X));
%         g(ix)=0;
%     end
function [f,g]=act2(X)
        f=0.1 +log(1+exp(X));
        ix= isinf(f);
        f(ix)=0.1;
        [m,im]=max(f,[],'all','linear');
        f(ix)=m;
        g=1./(1+exp(-X));
        g(ix)=g(im);
    end
%     function [f,g]=linear(W,X)
%         f=W*X;
%         g=nan([size(W,1),size(W,2),size(X,2)]);
%         for ix=1:size(W,1)
%             g(ix,:,:)=X;
%         end
%     end
    function ix = getIX(n)
        ix=1:n;
    end
    function [f,g]= gaussianLoss(X,y,p,mu)
        sigmaSq=X(1,:);
        %sigmaSq=X(2,:);
        L1= 0.5 * log (2*pi* sigmaSq) + 0.5 * ((y-mu).^2)./sigmaSq;
        f= mean(p .* L1);
        gSS= 0.5./sigmaSq - 0.5 * ((y-mu).^2)./(sigmaSq.^2);
        %gSS = mean(p.*gSS1);
        %gmu = - (y-mu)./sigmaSq;
        %g=[gmu;gSS];
        g=gSS;
        if any(isnan(g),'all')
            disp('hi');
        end
        %gmu = mean(p.*gmu1);
    end
end

