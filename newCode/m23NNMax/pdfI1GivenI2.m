function p = pdfI1GivenI2(I1,I2,zeta)
   diff=I1-I2;
   ix= diff>=0;
   p=zeros(length(I1),1);
   [mu,ss] = predictMuAndSS(zeta.D1,I2(ix));
   p(ix) = tnPdf(diff(ix),mu,sqrt(ss));
end

