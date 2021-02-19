function p = snPdfOL(x,mu,omega,lambda)
p = (2/omega) * normpdf((x-mu)/omega) .* normcdf(lambda * (x-mu)/omega);
end

