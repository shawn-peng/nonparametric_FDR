function [v,w] = VW(x,theta)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here\
mu=theta.mu;Gamma=theta.Gamma;Delta=theta.Delta;
[omega,~,delta]=GD2OL(Gamma,Delta);
loc=delta/omega*(x-mu);
scale=sqrt(1-delta^2);
[v,w]=momentsTN(loc,scale);
end

