function [f,g] = softPl(X)
        f=0.1 +log(1+exp(X));
        ix= isinf(f);
        f(ix)=0.1;
        [m,im]=max(f,[],'all','linear');
        f(ix)=m;
        g=1./(1+exp(-X));
        g(ix)=g(im);
end
