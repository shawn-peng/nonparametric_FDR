function [f,g] = relu(X)
     f=X;
     ix=X<0;  
     f(ix)= X(ix)*0.1;
    %f(ix)=0;
     g=ones(size(X));
     g(ix)=0.1;
end

