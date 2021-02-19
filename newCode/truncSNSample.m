function [smp,wt] = truncSNSample(theta,lb,ub)
n=10000;
xx = theta.mu + theta.Delta*abs(randn(n,1)) + sqrt(theta.Gamma)*randn(n,1);
for i = 1:n
    tree=makeTree(xx(i),tree,0);
end
for i=1:length(lb)
    p(i) = areaBetween(lb(i),ub(i));
end



%sigma=sigmaVec(thetaD2,z);

    function  tree = makeTree(x,tree,lSize)
        if isempty(tree)
            st=struct('x',x,'lSize',lSize+1);
            tree{1}=st;
            tree{2}={};
            tree{3}={};
        else
            st=tree{1};
            if (x < st.x)
                st.lSize=st.lSize+1;
                tree{2}=makeTree(x,tree{2},lSize);
            else
                %st.rSize=st.rSize+1;
                tree{3}=makeTree(x,tree{3},st.lSize);
            end
            tree{1}=st;
        end
    end
    function [a,s] = areaBetweenAndSample(l,u)
        if isfinite(l)
           stl=search(l,tree);
           c1=stl.lSize;
        else
           c1=0;
        end
        if isfinite(u)
           stu=search(u,tree);
           c2=stu.lSize;
        else
          c2=n;
        end
        a=(c2-c1)/n;
        
    end
end

 