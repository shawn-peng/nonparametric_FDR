function [s1Dens, s2Dens,S1,S2] = densEst_hist(s1,s2,~)
%Estimates the density of x and x1 using finite gaussian mixtures. The
%number of components is selected using AIC. The components obtained for x1
%are resused to fit x.
set(0,'DefaultFigureVisible', 'off');
if exist('histogram') == 2
    h1=histogram(s2);
    binWidth=h1.BinWidth;
else
    [~,centers]=hist(s2);
    binWidth=centers(2)-centers(1);
end
set(0,'DefaultFigureVisible', 'on');
%figure('visible','off');
xmin=min([s2;s1]);
xmax=max([s2;s1]);
binEdges=xmin:binWidth:(xmax+binWidth);
binEdges=binEdges(:);
numBins=length(binEdges)-1;
s2Dens=toUnifMixture(repmat(1/numBins,1,numBins), binEdges(1:numBins),binEdges(2:numBins+1)); 
S2=s2Dens.bin(s2);
s2Dens=s2Dens.fixedCompsFit(S2);
%s2Dens=fixedCompsFit(s2Dens,s2);
s1Dens=toUnifMixture(repmat(1/numBins,1,numBins), binEdges(1:numBins),binEdges(2:numBins+1)); 
S1=s1Dens.bin(s1);
s1Dens=s1Dens.fixedCompsFit(S1);
%s1Dens=fixedCompsFit(s1Dens,s1);
%x=sort(mixSample);
%x1=sort(comp1Sample);
%hold off
%plot(x,pdf(mixDist,x));
  %   hold on;     
%plot(x1,pdf(comp1Dist,x1));
%comp1Dist=fixedCompsFit(comp1Dist,x1);
%plot(x,pdf(comp1Dist,x));
     %hold off;
end