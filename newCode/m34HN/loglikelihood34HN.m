function ll = loglikelihood34HN(s1,s2,s3,zeta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[~,~,~,p]=pVec34HN(s1,s2,s3,zeta);
p(p==0)=min(p(p~=0));
ll=mean(log(p));
end

