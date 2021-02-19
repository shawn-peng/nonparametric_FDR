function theta = momEst(X)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
m1 = mean(X);
m2 = var(X);
m3 = skewness(X);
a1 = sqrt(2/pi);
b1 = (4/pi - 1) / a1;

delta = sign(m3)*min(1-10^-7,  1/ sqrt(a1^2 + m2 .* (b1 / abs(m3)).^(2/3)));

omega = sqrt(m2 ./ (1 - a1^2 * delta.^2));

mu = m1 - a1*delta*omega;

lambda = delta/sqrt(1-delta^2);
theta.mu=mu;
[theta.Gamma,theta.Delta]=OL2GD(omega,lambda);
theta.Gamma=max(theta.Gamma,1);
end

