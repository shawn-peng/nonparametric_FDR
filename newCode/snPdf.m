function p = snPdf(x,theta)
mu=theta.mu;Gamma=theta.Gamma;Delta=theta.Delta;
[omega,lambda]=GD2OL(Gamma,Delta);
p=snPdfOL(x,mu,omega,lambda);
end

