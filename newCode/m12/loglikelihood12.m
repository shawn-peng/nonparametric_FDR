function ll = loglikelihood12(s1,zeta);
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[~,~,p]=pVec12(s1,zeta);
p(p==0)=min(p(p~=0));
ll=mean(log(p));
end

