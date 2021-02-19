function [omega,lambda,delta] = GD2OL(Gamma,Delta)
omega =sqrt(Gamma + Delta^2);
lambda = sign(Delta)*sqrt(Delta^2/Gamma);
delta=lambda/sqrt(1+lambda^2);
end

