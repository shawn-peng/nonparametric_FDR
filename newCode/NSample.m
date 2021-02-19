function xx = NSample(theta,n)
% xx = theta.mu + theta.sigma*randn(n, 1);
xx = normrnd(theta.mu, theta.sigma, n, 1);
end

