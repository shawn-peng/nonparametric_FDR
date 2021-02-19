function [cdf,x] = sampleCDF(x)
x=sort(x);
cdf=(1:length(x))/length(x);
cdf=cdf(:);
end

