function pp = tnPdf(x,mu,sigma)
p=normpdf(x,mu,sigma);
p(x<0)=0;
q=1-normcdf(zeros(size(x)),mu,sigma);
pp=p./q;
pp(isinf(pp))=max(pp(~isinf(pp)));
end