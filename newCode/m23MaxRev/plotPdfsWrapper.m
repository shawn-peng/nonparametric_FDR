function  plotPdfsWrapper(s1,s2,zeta)
 [pdfM,pdfC,pdfI1]=s1Dens(s1,zeta);
 %[pdfM_2,pdfC_2,pdfI1_1, pdfI1_2]=s2Dens(s2,zeta);
 %plotPdfs(s1,s2,pdfM,pdfC,pdfI1,pdfM_2,pdfC_2,pdfI1_1,pdfI1_2)
 subplot(2,2,3)
 figure;
 cdfT=empiricalCdf(s1,s1);
 [~,cdfE]=FDR(s1,zeta);
 plot(s1,cdfT);
 hold on
 plot(s1,cdfE);
 delta=areaBetweenCDF(cdfT,cdfE,s1);
 title(['\deltacdf:',num2str(delta,4)])
 hold off;
 
end

