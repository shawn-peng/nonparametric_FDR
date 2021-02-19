function [f,g]=sigmoid(X)
        f=1./(1+exp(-X));
        g=f.*(1-f);
end