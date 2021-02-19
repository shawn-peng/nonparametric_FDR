function ll = loglikelihood23HN(s1,s2,zeta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[~,~,~,p]=pVec23HN(s1,s2,zeta);
p(p==0)=min(p(p~=0));
ll=mean(log(p));
end

