function pp = tnCdf(x,mu,sigma)
p0=normcdf(zeros(length(x),1),mu,sigma);
p=normcdf(x,mu,sigma);
pp=(p-p0)./(1-p0);
pp(isinf(pp))=max(pp(~isinf(pp)));
end
