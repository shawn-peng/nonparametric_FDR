function [Gamma,Delta,delta] = OL2GD(omega,lambda)
delta=lambda/sqrt(1+lambda^2);
Gamma = omega^2*(1-delta^2);
Delta = omega * delta;
end

