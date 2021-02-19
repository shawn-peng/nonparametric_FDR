function sigma = sigmaVec(theta,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%sigma=max(theta.a*x + theta.b,10^-7);
sigma=sim(theta.net,x');
sigma(sigma<=10^-7)=10^-7;
sigma=sigma';
end

