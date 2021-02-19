function [mu,ss] = predictMuAndSS(theta,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% ss=theta.netSS.predict(x);
% mu=theta.netMu.predict(x);
ss = max((5 + x * 10) .^ 2, 0.01);
mu = 0 * ones(size(x));
% ss= theta.ss + zeros(length(x),1);
% mu= theta.mu + zeros(length(x),1);
end

