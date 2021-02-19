function p = pdfI1GivenI2(zeta,I1,I2)
   diff=I1-I2;
   ix= diff>=0;
   I22=double(I2);
   p=zeros(length(I1),1);
   [mu,ss] = predictMuAndSS(zeta.D1,I22(ix));
   p(ix) = tnPdf(diff(ix),mu,sqrt(ss));
end

