function [mu,ss] = predictMuAndSS(theta,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%ss=theta.netSS.predict(x);
%mu=theta.netMu.predict(x);
ss= theta.ss + zeros(length(x),1);
mu= theta.mu + zeros(length(x),1);
end

