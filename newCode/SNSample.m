function xx = SNSample(theta,n)
xx = theta.mu + theta.Delta*abs(randn(n,1)) + sqrt(theta.Gamma)*randn(n,1);
end

