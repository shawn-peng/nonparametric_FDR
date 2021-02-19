function [mu,ss] = predictMuAndSS(theta,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ss=theta.netSS.predict(x);
mu=theta.netMu.predict(x);
end

